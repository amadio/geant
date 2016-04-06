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

  VECPHYS_CUDA_HEADER_HOST
  PhotoElectronSauterGavrila(Random_t* states, int threadId = -1);

  VECPHYS_CUDA_HEADER_BOTH
  PhotoElectronSauterGavrila(Random_t* states, int threadId,
                             GUAliasSampler* sampler);

  VECPHYS_CUDA_HEADER_BOTH
  ~PhotoElectronSauterGavrila() {}

  VECPHYS_CUDA_HEADER_HOST
  void Initialization();

  //interfaces for tables
  VECPHYS_CUDA_HEADER_HOST
  void BuildCrossSectionTablePerAtom(int Z);

  VECPHYS_CUDA_HEADER_HOST
  void BuildPdfTable(int Z, double *p);

  //Alternative Interact method to test energy dependent subtasks for a
  //specific model. Eeventually this method should replace the Interact
  //method of EmBaseModel

  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void ModelInteract(GUTrack&  projectile,
                     const int targetElement,
                     GUTrack&  secondary );

  //vector
#ifndef VECPHYS_NVCC
  template <typename Backend>
  void ModelInteract(GUTrack_v& inProjectile,
                     const int* targetElements,
                     GUTrack_v& outSecondaryV);
#endif


private:
  // Implementation methods
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  CrossSectionKernel(typename Backend::Double_t  energyIn,
                     typename Backend::Index_t   zElement);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void
  InteractKernel(typename Backend::Double_t energyIn,
                 typename Backend::Index_t  zElement,
                 typename Backend::Double_t& energyOut,
                 typename Backend::Double_t& sinTheta);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void
  InteractKernelCR(typename Backend::Double_t energyIn,
                   typename Backend::Index_t   zElement,
                   typename Backend::Double_t& energyOut,
                   typename Backend::Double_t& sinTheta);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void
  InteractKernelUnpack(typename Backend::Double_t energyIn,
                       typename Backend::Index_t   zElement,
                       typename Backend::Double_t& energyOut,
                       typename Backend::Double_t& sinTheta,
                       typename Backend::Bool_t &status);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  GetPhotoElectronEnergy(typename Backend::Double_t energyIn,
                         typename Backend::Index_t zElement);

  template<class Backend>
  inline
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  SampleSequential(typename Backend::Double_t A,
                   typename Backend::Double_t Ap2,
                   typename Backend::Double_t B,
                   typename Backend::Double_t grej) const;

  VECPHYS_CUDA_HEADER_BOTH
  void SampleByCompositionRejection(int    Z,
                                    double energyIn,
                                    double& energyOut,
                                    double& sinTheta);

  VECPHYS_CUDA_HEADER_BOTH double
  GetG4CrossSection(double  energyIn,
                    const int zElement);


  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSectionK( int Zelement, double Ein, double outEphoton ) const;

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<PhotoElectronSauterGavrila>;

  //private:

};

//Implementation

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
PhotoElectronSauterGavrila::CrossSectionKernel(typename Backend::Double_t energy,
                                               typename Backend::Index_t  Z)
{
  typedef typename Backend::Double_t Double_t;

  Double_t sigma = 0.;

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

  Double_t energy2 = energy*energy;
  Double_t energy3 = energy*energy2;
  Double_t energy4 = energy2*energy2;

  Double_t sgima = fSandiaCof[0]/energy  + fSandiaCof[1]/energy2 +
    fSandiaCof[2]/energy3 + fSandiaCof[3]/energy4;

  return sigma;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void PhotoElectronSauterGavrila::
InteractKernel(typename Backend::Double_t  energyIn,
               typename Backend::Index_t   zElement,
               typename Backend::Double_t& energyOut,
               typename Backend::Double_t& sinTheta)
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

  //energy of photo-electron: Sandia parameterization
  energyOut = GetPhotoElectronEnergy<Backend>(energyIn,zElement) ;

  //sample angular distribution of photo-electron

  Index_t   irow;
  Index_t   icol;
  Double_t  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,irow,icol,fraction);

  Double_t probNA;
  Index_t  aliasInd;

  Double_t ncol(fAliasSampler->GetSamplesPerEntry());
  Index_t   index = ncol*irow + icol;
  fAliasSampler->GatherAlias<Backend>(index,probNA,aliasInd);

  Double_t mininum = -1.0;
  Double_t deltaE = 2.0;

  Double_t cosTheta = mininum + fAliasSampler->SampleX<Backend>(deltaE,probNA,
                                                      aliasInd,icol,fraction);

  sinTheta = Sqrt((1+cosTheta)*(1-cosTheta));
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
PhotoElectronSauterGavrila::
GetPhotoElectronEnergy(typename Backend::Double_t energy,
                       typename Backend::Index_t  zElement)
{
  // this method is not vectorizable and only for the scalar backend

  typedef typename Backend::Int_t Int_t;
  typedef typename Backend::Double_t Double_t;

  // Photo electron energy
  Double_t energyOut = 0.;

  // Select atomic shell
  assert (zElement>0 && zElement <101);

  Int_t nShells = fNumberOfShells[zElement];

  Int_t i = 0;
  Double_t bindingEnergy =0;

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

#ifndef VECPHYS_NVCC
template<>
inline
VECPHYS_CUDA_HEADER_BOTH
typename kVc::Double_t
PhotoElectronSauterGavrila::
GetPhotoElectronEnergy<kVc>(typename kVc::Double_t energy,
                            typename kVc::Index_t  zElement)
{
  kVc::Double_t energyOut;

  for(int i = 0; i < kVc::kSize ; ++i) {
    energyOut[i] = GetPhotoElectronEnergy<kScalar>(energy[i],zElement[i]);
  }

  return energyOut;
}

#endif


template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void
PhotoElectronSauterGavrila::InteractKernelCR(typename Backend::Double_t  energyIn,
                                             typename Backend::Index_t   zElement,
                                             typename Backend::Double_t& energyOut,
                                             typename Backend::Double_t& sinTheta)
{
  typedef typename Backend::Double_t Double_t;
  //  typedef typename Backend::Bool_t Bool_t;

  //energy of photo-electron: Sandia parameterization
  energyOut = GetPhotoElectronEnergy<Backend>(energyIn,zElement) ;

  //sample angular direction according to SauterGavrilaAngularDistribution
  Double_t tau = energyIn/electron_mass_c2;

  /*
  const double taulimit = 50.0;
  Bool_t highE = tau > taulimit;
  cosTheta = 1.0;
  if(Backend::early_returns && IsFull(highE)) return;
  */

  Double_t gamma     = tau + 1.0;
  Double_t beta      = sqrt(tau*(tau + 2.0))/gamma;

  Double_t A = (1-beta)/beta;
  Double_t Ap2 = A + 2;
  Double_t B   = 0.5*beta*gamma*(gamma-1.)*(gamma-2.);
  Double_t grej = 2*(1+A*B)/A;

  Double_t z = SampleSequential<Backend>(A,Ap2,B,grej);

  //  MaskedAssign(!highE, 1.0 - z , &cosTheta);
  sinTheta = Sqrt(z*(2-z)); // cosTheta = 1 -z

}
template<class Backend>
inline
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
PhotoElectronSauterGavrila::SampleSequential(typename Backend::Double_t A,
                                             typename Backend::Double_t Ap2,
                                             typename Backend::Double_t B,
                                             typename Backend::Double_t grej) const
{
  typedef typename Backend::Double_t Double_t;

  Double_t z;
  Double_t g;

  do {
    Double_t q = UniformRandom<Backend>(fRandomState,fThreadId);
    z = 2*A*(2*q + Ap2*sqrt(q))/(Ap2*Ap2 - 4*q);
    g = (2 - z)*(1.0/(A + z) + B);
  } while(g < UniformRandom<Backend>(fRandomState,fThreadId)*grej);

  return z;
}

#ifndef VECPHYS_NVCC
template<>
inline
VECPHYS_CUDA_HEADER_BOTH
typename kVc::Double_t
PhotoElectronSauterGavrila::SampleSequential<kVc>(typename kVc::Double_t A,
                                                  typename kVc::Double_t Ap2,
                                                  typename kVc::Double_t B,
                                                  typename kVc::Double_t grej) const
{
  typedef typename kVc::Double_t Double_t;

  Double_t z;
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
VECPHYS_CUDA_HEADER_BOTH void
PhotoElectronSauterGavrila::InteractKernelUnpack(typename Backend::Double_t energyIn,
                                                 typename Backend::Index_t   zElement,
                                                 typename Backend::Double_t& energyOut,
                                                 typename Backend::Double_t& sinTheta,
                                                 typename Backend::Bool_t &status)
{
  //dummy for now
  energyOut = energyIn;
  sinTheta =  0;
}

//Alternative Interact method

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH
void PhotoElectronSauterGavrila::ModelInteract(GUTrack&  inProjectile,
                                               const int targetElement,
                                               GUTrack&  outSecondary )
{
  typedef typename Backend::Double_t Double_t;

  Double_t energyIn = inProjectile.E;

  //check for the validity of energy
  if(energyIn < fLowEnergyLimit || energyIn > fHighEnergyLimit) return;

  Double_t energyOut =0;
  Double_t sinTheta = 0;

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
#ifndef VECPHYS_NVCC

template <typename Backend>
void PhotoElectronSauterGavrila::ModelInteract(GUTrack_v& inProjectile,
                                               const int* targetElements,
                                               GUTrack_v& outSecondary)
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

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
  int numChunks= (nTracks/Double_t::Size);

  for(int i= 0; i < numChunks ; ++i) {
    Double_t energyIn(&inProjectile.E[ibase]);
    Double_t sinTheta(0.);
    Double_t energyOut;

    Index_t  zElement(targetElements[ibase]);

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

    ibase+= Double_t::Size;
  }

  //leftover - do scalar (temporary)
  for(int i = numChunks*Double_t::Size ; i < nTracks ; ++i) {

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
