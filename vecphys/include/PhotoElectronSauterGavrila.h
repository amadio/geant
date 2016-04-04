#ifndef PhotoElectronSauterGavrila_H
#define PhotoElectronSauterGavrila_H 1

#include "base/Global.h"
#include "base/PhysicalConstants.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "StaticSandiaData.h"

#include "EmModelBase.h"

#include <stdlib.h>

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class PhotoElectronSauterGavrila : public EmModelBase<PhotoElectronSauterGavrila>
{
public:

  VECCORE_CUDA_HOST
  PhotoElectronSauterGavrila(Random_t* states, int threadId = -1);

  VECCORE_CUDA_HOST_DEVICE
  PhotoElectronSauterGavrila(Random_t* states, int threadId,
                             GUAliasSampler* sampler);

  VECCORE_CUDA_HOST_DEVICE
  ~PhotoElectronSauterGavrila() {}

  VECCORE_CUDA_HOST
  void Initialization();

  //interfaces for tables
  VECCORE_CUDA_HOST
  void BuildCrossSectionTablePerAtom(int Z);

  VECCORE_CUDA_HOST
  void BuildPdfTable(int Z, double *p);

  //Alternative Interact method to test energy dependent subtasks for a
  //specific model. Eeventually this method should replace the Interact
  //method of EmBaseModel

  template <typename Backend>
  VECCORE_CUDA_HOST_DEVICE
  void ModelInteract(GUTrack&  projectile,
                     const int targetElement,
                     GUTrack&  secondary );

  //vector
#ifndef VECCORE_NVCC
  template <typename Backend>
  void ModelInteract(GUTrack_v& inProjectile,
                     const int* targetElements,
                     GUTrack_v& outSecondaryV);
#endif


private:
  // Implementation methods
  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  CrossSectionKernel(typename Backend::Double_v  energyIn,
                     Index_v<typename Backend::Double_v>   zElement);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE void
  InteractKernel(typename Backend::Double_v energyIn,
                 Index_v<typename Backend::Double_v>  zElement,
                 typename Backend::Double_v& energyOut,
                 typename Backend::Double_v& sinTheta);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE void
  InteractKernelCR(typename Backend::Double_v energyIn,
                   Index_v<typename Backend::Double_v>   zElement,
                   typename Backend::Double_v& energyOut,
                   typename Backend::Double_v& sinTheta);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE void
  InteractKernelUnpack(typename Backend::Double_v energyIn,
                       Index_v<typename Backend::Double_v>   zElement,
                       typename Backend::Double_v& energyOut,
                       typename Backend::Double_v& sinTheta,
                       Mask_v<typename Backend::Double_v> &status);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  GetPhotoElectronEnergy(typename Backend::Double_v energyIn,
                         Index_v<typename Backend::Double_v> zElement);

  template<class Backend>
  inline
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  SampleSequential(typename Backend::Double_v A,
                   typename Backend::Double_v Ap2,
                   typename Backend::Double_v B,
                   typename Backend::Double_v grej) const;

  VECCORE_CUDA_HOST_DEVICE
  void SampleByCompositionRejection(int    Z,
                                    double energyIn,
                                    double& energyOut,
                                    double& sinTheta);

  VECCORE_CUDA_HOST_DEVICE double
  GetG4CrossSection(double  energyIn,
                    const int zElement);


  VECCORE_CUDA_HOST_DEVICE
  double CalculateDiffCrossSectionK( int Zelement, double Ein, double outEphoton ) const;

  VECCORE_CUDA_HOST_DEVICE
  double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<PhotoElectronSauterGavrila>;

  //private:

};

//Implementation

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
PhotoElectronSauterGavrila::CrossSectionKernel(typename Backend::Double_v energy,
                                               Index_v<typename Backend::Double_v>  Z)
{
  using Double_v = typename Backend::Double_v;

  Double_v sigma = 0.;

  //Sandia parameterization for Z < 100
  //  int Z = zElement;

  int    fCumulInterval[101]  = {0};
  double fSandiaCof[4]        = {0.0};

  fCumulInterval[0] = 1;

  //scan - move it to constructor or use pre-built table
  for (int iz = 1; iz < 101; ++iz) {
    fCumulInterval[iz] = fCumulInterval[iz-1] + fNbOfIntervals[iz];
  }

  double Emin  = fSandiaTable[fCumulInterval[Z-1]][0]*keV;

  int interval = fNbOfIntervals[Z] - 1;
  int row = fCumulInterval[Z-1] + interval;

  while ((interval>0) && (energy<fSandiaTable[row][0]*keV)) {
    --interval;
    row = fCumulInterval[Z-1] + interval;
  }

  if (energy >= Emin) {
    double AoverAvo = Z*amu/fZtoAratio[Z];
    fSandiaCof[0]=AoverAvo*funitc[1]*fSandiaTable[row][1];
    fSandiaCof[1]=AoverAvo*funitc[2]*fSandiaTable[row][2];
    fSandiaCof[2]=AoverAvo*funitc[3]*fSandiaTable[row][3];
    fSandiaCof[3]=AoverAvo*funitc[4]*fSandiaTable[row][4];
  }
  else {
    fSandiaCof[0] = fSandiaCof[1] = fSandiaCof[2] = fSandiaCof[3] = 0.;
  }

  Double_v energy2 = energy*energy;
  Double_v energy3 = energy*energy2;
  Double_v energy4 = energy2*energy2;

  Double_v sgima = fSandiaCof[0]/energy  + fSandiaCof[1]/energy2 +
    fSandiaCof[2]/energy3 + fSandiaCof[3]/energy4;

  return sigma;
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
void PhotoElectronSauterGavrila::
InteractKernel(typename Backend::Double_v  energyIn,
               Index_v<typename Backend::Double_v>   zElement,
               typename Backend::Double_v& energyOut,
               typename Backend::Double_v& sinTheta)
{
  typedef Index_v<typename Backend::Double_v>  Index_v<Double_v>;
  using Double_v = typename Backend::Double_v;

  //energy of photo-electron: Sandia parameterization
  energyOut = GetPhotoElectronEnergy<Backend>(energyIn,zElement) ;

  //sample angular distribution of photo-electron

  Index_v<Double_v>   irow;
  Index_v<Double_v>   icol;
  Double_v  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,irow,icol,fraction);

  Double_v probNA;
  Index_v<Double_v>  aliasInd;

  Double_v ncol(fAliasSampler->GetSamplesPerEntry());
  Index_v<Double_v>   index = ncol*irow + icol;
  fAliasSampler->GatherAlias<Backend>(index,probNA,aliasInd);

  Double_v mininum = -1.0;
  Double_v deltaE = 2.0;

  Double_v cosTheta = mininum + fAliasSampler->SampleX<Backend>(deltaE,probNA,
                                                      aliasInd,icol,fraction);

  sinTheta = Sqrt((1+cosTheta)*(1-cosTheta));
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
PhotoElectronSauterGavrila::
GetPhotoElectronEnergy(typename Backend::Double_v energy,
                       Index_v<typename Backend::Double_v>  zElement)
{
  // this method is not vectorizable and only for the scalar backend

  typedef typename Backend::Int_t Int_t;
  using Double_v = typename Backend::Double_v;

  // Photo electron energy
  Double_v energyOut = 0.;

  // Select atomic shell
  assert (zElement>0 && zElement <101);

  Int_t nShells = fNumberOfShells[zElement];

  Int_t i = 0;
  Double_v bindingEnergy =0;

  for( ; i < nShells ; ++i) {
    bindingEnergy = fBindingEnergies[fIndexOfShells[zElement] + i]*eV;
    if(energy >= bindingEnergy ) { break; }
  }

  // Normally one shell is available
  if (i < nShells) {
    bindingEnergy = fBindingEnergies[fIndexOfShells[zElement] + i]*eV;

    // update by deexcitation goes here

    energyOut = energy - bindingEnergy;
  }

  return energyOut;
}

#ifndef VECCORE_NVCC
template<>
inline
VECCORE_CUDA_HOST_DEVICE
typename kVc::Double_v
PhotoElectronSauterGavrila::
GetPhotoElectronEnergy<kVc>(typename kVc::Double_v energy,
                            typename kVc::Index_v<Double_v>  zElement)
{
  kVc::Double_v energyOut;

  for(int i = 0; i < kVc::kSize ; ++i) {
    energyOut[i] = GetPhotoElectronEnergy<kScalar>(energy[i],zElement[i]);
  }

  return energyOut;
}

#endif


template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
PhotoElectronSauterGavrila::InteractKernelCR(typename Backend::Double_v  energyIn,
                                             Index_v<typename Backend::Double_v>   zElement,
                                             typename Backend::Double_v& energyOut,
                                             typename Backend::Double_v& sinTheta)
{
  using Double_v = typename Backend::Double_v;
  //  typedef Mask_v<typename Backend::Double_v> Mask_v<Double_v>;

  //energy of photo-electron: Sandia parameterization
  energyOut = GetPhotoElectronEnergy<Backend>(energyIn,zElement) ;

  //sample angular direction according to SauterGavrilaAngularDistribution
  Double_v tau = energyIn/electron_mass_c2;

  /*
  const double taulimit = 50.0;
  Mask_v<Double_v> highE = tau > taulimit;
  cosTheta = 1.0;
  if(Backend::early_returns && IsFull(highE)) return;
  */

  Double_v gamma     = tau + 1.0;
  Double_v beta      = sqrt(tau*(tau + 2.0))/gamma;

  Double_v A = (1-beta)/beta;
  Double_v Ap2 = A + 2;
  Double_v B   = 0.5*beta*gamma*(gamma-1.)*(gamma-2.);
  Double_v grej = 2*(1+A*B)/A;

  Double_v z = SampleSequential<Backend>(A,Ap2,B,grej);

  //  MaskedAssign(!highE, 1.0 - z , &cosTheta);
  sinTheta = Sqrt(z*(2-z)); // cosTheta = 1 -z

}
template<class Backend>
inline
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
PhotoElectronSauterGavrila::SampleSequential(typename Backend::Double_v A,
                                             typename Backend::Double_v Ap2,
                                             typename Backend::Double_v B,
                                             typename Backend::Double_v grej) const
{
  using Double_v = typename Backend::Double_v;

  Double_v z;
  Double_v g;

  do {
    Double_v q = UniformRandom<Backend>(fRandomState,fThreadId);
    z = 2*A*(2*q + Ap2*sqrt(q))/(Ap2*Ap2 - 4*q);
    g = (2 - z)*(1.0/(A + z) + B);
  } while(g < UniformRandom<Backend>(fRandomState,fThreadId)*grej);

  return z;
}

#ifndef VECCORE_NVCC
template<>
inline
VECCORE_CUDA_HOST_DEVICE
typename kVc::Double_v
PhotoElectronSauterGavrila::SampleSequential<kVc>(typename kVc::Double_v A,
                                                  typename kVc::Double_v Ap2,
                                                  typename kVc::Double_v B,
                                                  typename kVc::Double_v grej) const
{
  typedef typename kVc::Double_v Double_v;

  Double_v z;
  double g;

  for(int i = 0; i < kVc::kSize ; ++i) {
    do {
      double q = UniformRandom<kScalar>(fRandomState,fThreadId);
      z[i] = 2*A[i]*(2*q + Ap2[i]*sqrt(q))/(Ap2[i]*Ap2[i] - 4*q);
      g = (2 - z[i])*(1.0/(A[i] + z[i]) + B[i]);
    } while(g < UniformRandom<kScalar>(fRandomState,fThreadId)*grej[i]);
  }

  return z;
}
#endif

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
PhotoElectronSauterGavrila::InteractKernelUnpack(typename Backend::Double_v energyIn,
                                                 Index_v<typename Backend::Double_v>   zElement,
                                                 typename Backend::Double_v& energyOut,
                                                 typename Backend::Double_v& sinTheta,
                                                 Mask_v<typename Backend::Double_v> &status)
{
  //dummy for now
  energyOut = energyIn;
  sinTheta =  0;
}

//Alternative Interact method

template <typename Backend>
VECCORE_CUDA_HOST_DEVICE
void PhotoElectronSauterGavrila::ModelInteract(GUTrack&  inProjectile,
                                               const int targetElement,
                                               GUTrack&  outSecondary )
{
  using Double_v = typename Backend::Double_v;

  Double_v energyIn = inProjectile.E;

  //check for the validity of energy
  if(energyIn < fLowEnergyLimit || energyIn > fHighEnergyLimit) return;

  Double_v energyOut =0;
  Double_v sinTheta = 0;

  //a good upper bound of photon energy to apply the alias method for
  //the SauterGavrila angular distribution (above this energy region,
  //dsigma/dcos(theta) is squeezed toward 1
  const double aliaslimit = 1.0*MeV;

  //lower bound for the approximation dsigma/dcos(theta) =1 driven by Geant4
  //(note that the (geant4) composition and rejection method is very inefficient
  //for the model above this region)
  const double taulimit = 50.*electron_mass_c2;

  if(energyIn < aliaslimit) {
    InteractKernel<Backend>(energyIn,targetElement,energyOut,sinTheta);
  }
  else if(energyIn < taulimit) {
    InteractKernelCR<Backend>(energyIn,targetElement,energyOut,sinTheta);
  }
  else {
    energyOut = GetPhotoElectronEnergy<Backend>(energyIn,targetElement);
    sinTheta = 0; //cosTheta = 1.0;
  }

  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
}
#ifndef VECCORE_NVCC

template <typename Backend>
void PhotoElectronSauterGavrila::ModelInteract(GUTrack_v& inProjectile,
                                               const int* targetElements,
                                               GUTrack_v& outSecondary)
{
  typedef Index_v<typename Backend::Double_v>  Index_v<Double_v>;
  using Double_v = typename Backend::Double_v;

  //filtering energy regions for sampling methods - setable if necessary
  const double aliaslimit = 1.0*MeV;
  const double taulimit = 50.*electron_mass_c2;

  int nTracks = inProjectile.numTracks;
  double* start = inProjectile.E;
  auto indexAliasLimit = std::lower_bound(start,start+nTracks,aliaslimit) - start;
  auto indexTauLimit = std::lower_bound(start,start+nTracks,taulimit) - start;

  for(int j = 0; j < nTracks  ; ++j) {
    assert( (targetElements[j] > 0)  && (targetElements[j] <= maximumZ ) );
  }

  int ibase= 0;
  int numChunks= (nTracks/Double_v::Size);

  for(int i= 0; i < numChunks ; ++i) {
    Double_v energyIn(&inProjectile.E[ibase]);
    Double_v sinTheta(0.);
    Double_v energyOut;

    Index_v<Double_v>  zElement(targetElements[ibase]);

    if(ibase < indexAliasLimit) {
      InteractKernel<Backend>(energyIn,zElement,energyOut,sinTheta);
    }
    else if (ibase < indexTauLimit) {
      InteractKernelCR<Backend>(energyIn,zElement,energyOut,sinTheta);
    }
    else {
      energyOut = GetPhotoElectronEnergy<Backend>(energyIn,zElement);
      sinTheta = 0.; //cosTheta = 1.0;
    }

    ConvertXtoFinalState<Backend>(energyIn, energyOut, sinTheta, ibase, inProjectile, outSecondary);

    ibase+= Double_v::Size;
  }

  //leftover - do scalar (temporary)
  for(int i = numChunks*Double_v::Size ; i < nTracks ; ++i) {

    double senergyIn= inProjectile.E[i];
    double senergyOut, ssinTheta;
    //use InteractKernel for any leftover to be consistent with EmBaseModel
    InteractKernel<kScalar>(senergyIn,targetElements[i],senergyOut,ssinTheta);
    ConvertXtoFinalState_Scalar<kScalar>(senergyIn, senergyOut, ssinTheta, i, inProjectile, outSecondary);
  }
}

#endif


} // end namespace impl
} // end namespace vecphys

#endif
