#ifndef IonisationMoller_H
#define IonisationMoller_H 1

#include "base/Global.h"
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

  VECPHYS_CUDA_HEADER_HOST
  void Initialization();

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
  inline
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  SampleSequential(typename Backend::Double_t xmin,
                   typename Backend::Double_t xmax,
                   typename Backend::Double_t gg) const;

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
  double fDeltaRayThreshold; //delta-ray threshold
};

//Implementation
template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
IonisationMoller::CrossSectionKernel(typename Backend::Double_t energy,
                                     typename Backend::Index_t  Z)
{
  //the total cross section for Moller scattering per atom
  //energy = kinetic energy

  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  Double_t sigmaOut = 0.;
  Bool_t belowLimit = Bool_t(false);
  //low energy limit
  belowLimit |= ( energy < fLowEnergyLimit );
  if(Backend::early_returns && IsFull(belowLimit)) return;

  //delta-ray cuff-off (material dependent) use 1.0*keV temporarily
  Double_t tmin = fDeltaRayThreshold;
  Double_t tmax = 0.5*energy;

  Double_t xmin  = tmin/energy;
  Double_t xmax  = tmax/energy;
  Double_t tau   = energy/electron_mass_c2;
  Double_t gam   = tau + 1.0;
  Double_t gamma2= gam*gam;
  Double_t beta2 = tau*(tau + 2)/gamma2;

  //Moller (e-e-) scattering
  //H. Messel and D.F. Crawford, Pergamon Press, Oxford (1970)
  //G4MollerBhabhaModel::ComputeCrossSectionPerAtom

  Double_t gg = (2.0*gam - 1.0)/gamma2;
  sigmaOut = ((xmax - xmin)*(1.0 - gg + 1.0/(xmin*xmax)
          + 1.0/((1.0-xmin)*(1.0 - xmax)))
          - gg*Log( xmax*(1.0 - xmin)/(xmin*(1.0 - xmax)) ) ) / beta2;

  sigmaOut *= Z*twopi_mc2_rcl2/energy;

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
  Index_t  aliasInd;

  //this did not used to work - Fixed SW
  Double_t ncol(fAliasSampler->GetSamplesPerEntry());
  Index_t   index = ncol*irow + icol;
  fAliasSampler->GatherAlias<Backend>(index,probNA,aliasInd);

  Double_t mininumE = fDeltaRayThreshold;
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
  Double_t totalMomentum = Sqrt(energyIn*(energyIn + 2.0*electron_mass_c2));

  Double_t deltaMomentum = Sqrt(energyOut * (energyOut + 2.0*electron_mass_c2));
  Double_t cost =  energyOut * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);

  Double_t sint2 = (1.0 - cost)*(1. + cost);
  Double_t sinTheta;
  Bool_t condition2 = sint2 < 0.0;
  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );

  return sinTheta;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void
IonisationMoller::InteractKernelCR(typename Backend::Double_t  kineticEnergy,
                                   typename Backend::Index_t   zElement,
                                   typename Backend::Double_t& deltaKinEnergy,
                                   typename Backend::Double_t& sinTheta)
{
  typedef typename Backend::Bool_t Bool_t;
  typedef typename Backend::Double_t Double_t;

  //temporary - set by material
  Double_t cutEnergy = fDeltaRayThreshold;
  Double_t maxEnergy = 1.0*TeV;

  //based on G4MollerBhabhaModel::SampleSecondaries
  Double_t tmin = cutEnergy;
  Double_t tmax = 0.5*kineticEnergy;

  Bool_t condCut = (tmax < maxEnergy);
  MaskedAssign(!condCut, maxEnergy, &tmax);

  condCut |= (tmax >= tmin );

  if(Backend::early_returns && IsEmpty(condCut)) return;

  Double_t energy = kineticEnergy + electron_mass_c2;
  Double_t xmin   = tmin/kineticEnergy;
  Double_t xmax   = tmax/kineticEnergy;
  Double_t gam    = energy/electron_mass_c2;
  Double_t gamma2 = gam*gam;

  //Moller (e-e-) scattering
  Double_t gg = (2.0*gam - 1.0)/gamma2;

  Double_t x = SampleSequential<Backend>(xmin,xmax,gg);

  deltaKinEnergy = x * kineticEnergy;

  Double_t totalMomentum = Sqrt(kineticEnergy*(kineticEnergy + 2.0*electron_mass_c2));

  Double_t deltaMomentum =
    Sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  Double_t cost = deltaKinEnergy * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum );

  Bool_t condCos = (cost <= 1.0);
  MaskedAssign(!condCos, 1.0, &cost);

  Double_t sint2 = (1.0 - cost)*(1.0 + cost);

  Bool_t condSin2 = (sint2 >= 0.0);
  Double_t zero(0.0);
  CondAssign(condSin2, Sqrt(sint2), zero, &sinTheta);
}

template<class Backend>
inline
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
IonisationMoller::SampleSequential(typename Backend::Double_t xmin,
                                   typename Backend::Double_t xmax,
                                   typename Backend::Double_t gg) const
{
  typedef typename Backend::Int_t Int_t;
  typedef typename Backend::Double_t Double_t;

  Double_t  q;
  Double_t  x;
  Double_t  z;

  Double_t  y = 1.0 - xmax;
  Double_t grej =  1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));
  do {
    q = UniformRandom<Backend>(fRandomState,Int_t(fThreadId));
    x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
    y = 1.0 - x;
    z = 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y));
  } while(grej * UniformRandom<Backend>(fRandomState,Int_t(fThreadId)) > z);

  return x;
}

#ifndef VECPHYS_NVCC
template<>
inline
VECPHYS_CUDA_HEADER_BOTH
typename kVc::Double_t
IonisationMoller::SampleSequential<kVc>(typename kVc::Double_t xmin,
                                        typename kVc::Double_t xmax,
                                        typename kVc::Double_t gg) const
{
  typedef typename kVc::Double_t Double_t;

  Double_t  x;
  Double_t  y = 1.0 - xmax;
  Double_t  grej =  1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));

  double  q;
  double  z;

  for(int i = 0; i < kVc::kSize ; ++i) {
    do {
      q = UniformRandom<kScalar>(fRandomState,fThreadId);
      x[i] = xmin[i]*xmax[i]/(xmin[i]*(1.0 - q) + xmax[i]*q);
      y[i] = 1.0 - x[i];
      z = 1.0 - gg[i]*x[i] + x[i]*x[i]*(1.0 - gg[i] + (1.0 - gg[i]*y[i])/(y[i]*y[i]));
    } while(grej[i] * UniformRandom<kScalar>(fRandomState,fThreadId) > z);
  }
  return x;
}
#endif

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void
IonisationMoller::InteractKernelUnpack(typename Backend::Double_t energyIn,
                                       typename Backend::Index_t   zElement,
                                       typename Backend::Double_t& energyOut,
                                       typename Backend::Double_t& sinTheta,
                                       typename Backend::Bool_t &status)
{
  //dummy for now
  energyOut = energyIn;
  sinTheta =  0;
}

} // end namespace impl
} // end namespace vecphys

#endif
