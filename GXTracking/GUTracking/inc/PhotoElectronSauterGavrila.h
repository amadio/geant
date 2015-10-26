#ifndef PhotoElectronSauterGavrila_H
#define PhotoElectronSauterGavrila_H 1

#include "backend/Backend.h"
#include "base/PhysicalConstants.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "StaticSandiaData.h"

#include "EmModelBase.h"

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

  //interfaces for tables
  VECPHYS_CUDA_HEADER_HOST 
  void BuildCrossSectionTablePerAtom(int Z);

  VECPHYS_CUDA_HEADER_HOST
  void BuildPdfTable(int Z, double *p);

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
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t 
  GetPhotoElectronEnergy(typename Backend::Double_t energyIn,
                         typename Backend::Index_t zElement);

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
               typename Backend::Double_t& cosTheta)
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
  Double_t aliasInd;

  Index_t   index = fNcol*irow + icol;
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);
  
  Double_t mininum = -1.0;
  Double_t deltaE = 2.0;

  cosTheta = mininum + fAliasSampler->SampleX<Backend>(deltaE,probNA,
					     aliasInd,icol,fraction);

  //  sinTheta = sqrt((1+costTheta)*(1-costTheta));
  // @@@syj - rotation should be properly taken into account in the KernelInteract later
}    

/*
template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
void PhotoElectronSauterGavrila::
SampleByCompositionRejection(typename Backend::Double_t  energyIn,
                             typename Backend::Double_t& energyOut,
                             typename Backend::Double_t& cosTheta)
{
  typedef typename Backend::Double_t Double_t;

  //use the scalar implementation which is equivalent to Geant4
  energyOut = GetPhotoElectronEnergy<kScalar>(energyIn,10);

  //sample angular direction according to SauterGavrilaAngularDistribution

  Double_t tau = energyIn/electron_mass_c2;

  Double_t gamma     = tau + 1.0;
  Double_t beta      = sqrt(tau*(tau + 2.0))/gamma;

  Double_t A = (1-beta)/beta;
  Double_t Ap2 = A + 2;
  Double_t B   = 0.5*beta*gamma*(gamma-1.)*(gamma-2.);
  Double_t grej = 2*(1+A*B)/A;
  
  Double_t z;
  Double_t g;

  do { 
    Double_t q = UniformRandom<Backend>(fRandomState,fThreadId);
    z = 2*A*(2*q + Ap2*sqrt(q))/(Ap2*Ap2 - 4*q);
    g = (2 - z)*(1.0/(A + z) + B);
  } while(g < UniformRandom<Backend>(fRandomState,fThreadId)*grej);
  
  cosTheta = 1 - z;
}
*/

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

} // end namespace impl
} // end namespace vecphys

#endif
