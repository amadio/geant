#ifndef IonisationMoller_H
#define IonisationMoller_H 1

#include "backend/Backend.h"
#include "base/PhysicalConstants.h"

#include "GUConstants.h"
#include "GUTrack.h"

#include "EmModelBase.h"

namespace vecphys {

VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasSampler; )
inline namespace VECPHYS_IMPL_NAMESPACE {


class IonisationMoller : public EmModelBase<IonisationMoller>
{
public:

  VECPHYS_CUDA_HEADER_HOST
  IonisationMoller(Random_t* states = 0, int threadId = -1);

  VECPHYS_CUDA_HEADER_BOTH
  IonisationMoller(Random_t* states, int threadId, GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH
  ~IonisationMoller(){}

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
                 typename Backend::Index_t   zElement,
                 typename Backend::Double_t& energyOut,
                 typename Backend::Double_t& sinTheta);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  SampleSinTheta(typename Backend::Double_t energyIn,
                 typename Backend::Double_t energyOut) const; 

  VECPHYS_CUDA_HEADER_BOTH 
  void SampleByCompositionRejection(int    Z,
                                    double energyIn,
                                    double& energyOut,
                                    double& sinTheta);

  VECPHYS_CUDA_HEADER_BOTH double
  GetG4CrossSection(double  energyIn, 
                    const int zElement);
  
  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<IonisationMoller>;

private:
};

//Implementation
template<class Backend>
VECPHYS_CUDA_HEADER_BOTH 
typename Backend::Double_t
IonisationMoller::CrossSectionKernel(typename Backend::Double_t energy, 
                                     typename Backend::Index_t  Z)
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  Double_t sigmaOut = 0.;
  Bool_t belowLimit = Bool_t(false);
  //low energy limit
  belowLimit |= ( energy < fLowEnergyLimit );
  if(Backend::early_returns && IsFull(belowLimit)) return;  

  //delta-ray cuff-off (material dependent) use 1.0*keV temporarily
  Double_t tmin = 1.0*keV;
  Double_t tmax = 0.5*energy;
  
  Double_t xmin  = tmin/energy;
  Double_t xmax  = tmax/energy;
  Double_t tau   = energy/electron_mass_c2;
  Double_t gam   = tau + 1.0;
  Double_t gamma2= gam*gam;
  Double_t beta2 = tau*(tau + 2)/gamma2;
 
  //Moller (e-e-) scattering
 
  Double_t gg = (2.0*gam - 1.0)/gamma2;
  sigmaOut = ((xmax - xmin)*(1.0 - gg + 1.0/(xmin*xmax)
          + 1.0/((1.0-xmin)*(1.0 - xmax)))
          - gg*Log( xmax*(1.0 - xmin)/(xmin*(1.0 - xmax)) ) ) / beta2;
  
  //this is the case if one of E < belowLimit 
  MaskedAssign(belowLimit, 0.0,&sigmaOut);
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
IonisationMoller::InteractKernel(typename Backend::Double_t  energyIn, 
                                 typename Backend::Index_t   zElement,
                                 typename Backend::Double_t& energyOut,
                                 typename Backend::Double_t& sinTheta)
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

  Index_t   irow;
  Index_t   icol;
  Double_t  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,irow,icol,fraction);

  Double_t probNA;
  Double_t aliasInd;

  //this did not used to work - Fixed SW
  Index_t   index = fNcol*irow + icol;
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);
  
  Double_t mininumE = 0.1*keV;
  Double_t deltaE = energyIn/2.0 - mininumE;

  energyOut = mininumE 
    + fAliasSampler->SampleX<Backend>(deltaE,probNA,aliasInd,icol,fraction);
  sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t 
IonisationMoller::SampleSinTheta(typename Backend::Double_t energyIn,
                                 typename Backend::Double_t energyOut) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  //angle of the scatterred electron

  Double_t energy = energyIn + electron_mass_c2;
  Double_t totalMomentum = sqrt(energyIn*(energy + electron_mass_c2));

  Double_t deltaMomentum = sqrt(energyOut * (energyOut + 2.0*electron_mass_c2));
  Double_t cost =  energyOut * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);
  Double_t sint2 = (1.0 - cost)*(1. + cost);

  Double_t sinTheta = 0.5;
  Bool_t condition2 = sint2 < 0.0;

  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );   
  
  return sinTheta;
}

} // end namespace impl
} // end namespace vecphys

#endif
