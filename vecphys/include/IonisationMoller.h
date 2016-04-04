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

  VECCORE_CUDA_HOST
  IonisationMoller(Random_t* states = 0, int threadId = -1);

  VECCORE_CUDA_HOST_DEVICE
  IonisationMoller(Random_t* states, int threadId, GUAliasSampler* sampler);

  VECCORE_CUDA_HOST_DEVICE
  ~IonisationMoller(){}

  VECCORE_CUDA_HOST
  void Initialization();

  //interfaces for tables
  VECCORE_CUDA_HOST
  void BuildCrossSectionTablePerAtom(int Z);

  VECCORE_CUDA_HOST
  void BuildPdfTable(int Z, double *p);

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
                 Index_v<typename Backend::Double_v>   zElement,
                 typename Backend::Double_v& energyOut,
                 typename Backend::Double_v& sinTheta);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  SampleSinTheta(typename Backend::Double_v energyIn,
                 typename Backend::Double_v energyOut) const;

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
  inline
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  SampleSequential(typename Backend::Double_v xmin,
                   typename Backend::Double_v xmax,
                   typename Backend::Double_v gg) const;

  VECCORE_CUDA_HOST_DEVICE
  void SampleByCompositionRejection(int    Z,
                                    double energyIn,
                                    double& energyOut,
                                    double& sinTheta);

  VECCORE_CUDA_HOST_DEVICE double
  GetG4CrossSection(double  energyIn,
                    const int zElement);

  VECCORE_CUDA_HOST_DEVICE
  double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<IonisationMoller>;

private:
  double fDeltaRayThreshold; //delta-ray threshold
};

//Implementation
template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
IonisationMoller::CrossSectionKernel(typename Backend::Double_v energy,
                                     Index_v<typename Backend::Double_v>  Z)
{
  //the total cross section for Moller scattering per atom
  //energy = kinetic energy

  typedef Mask_v<typename Backend::Double_v>   Bool_t;
  using Double_v = typename Backend::Double_v;

  Double_v sigmaOut = 0.;
  Bool_t belowLimit = Bool_t(false);
  //low energy limit
  belowLimit |= ( energy < fLowEnergyLimit );
  if(Backend::early_returns && IsFull(belowLimit)) return;

  //delta-ray cuff-off (material dependent) use 1.0*keV temporarily
  Double_v tmin = fDeltaRayThreshold;
  Double_v tmax = 0.5*energy;

  Double_v xmin  = tmin/energy;
  Double_v xmax  = tmax/energy;
  Double_v tau   = energy/electron_mass_c2;
  Double_v gam   = tau + 1.0;
  Double_v gamma2= gam*gam;
  Double_v beta2 = tau*(tau + 2)/gamma2;

  //Moller (e-e-) scattering
  //H. Messel and D.F. Crawford, Pergamon Press, Oxford (1970)
  //G4MollerBhabhaModel::ComputeCrossSectionPerAtom

  Double_v gg = (2.0*gam - 1.0)/gamma2;
  sigmaOut = ((xmax - xmin)*(1.0 - gg + 1.0/(xmin*xmax)
          + 1.0/((1.0-xmin)*(1.0 - xmax)))
          - gg*Log( xmax*(1.0 - xmin)/(xmin*(1.0 - xmax)) ) ) / beta2;

  sigmaOut *= Z*twopi_mc2_rcl2/energy;

  //this is the case if one of E < belowLimit
  MaskedAssign(belowLimit, 0.0,&sigmaOut);
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
IonisationMoller::InteractKernel(typename Backend::Double_v  energyIn,
                                 Index_v<typename Backend::Double_v>   zElement,
                                 typename Backend::Double_v& energyOut,
                                 typename Backend::Double_v& sinTheta)
{
  typedef Index_v<typename Backend::Double_v>  Index_t;
  using Double_v = typename Backend::Double_v;

  Index_t   irow;
  Index_t   icol;
  Double_v  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,irow,icol,fraction);

  Double_v probNA;
  Index_t  aliasInd;

  //this did not used to work - Fixed SW
  Double_v ncol(fAliasSampler->GetSamplesPerEntry());
  Index_t   index = ncol*irow + icol;
  fAliasSampler->GatherAlias<Backend>(index,probNA,aliasInd);

  Double_v mininumE = fDeltaRayThreshold;
  Double_v deltaE = energyIn/2.0 - mininumE;

  energyOut = mininumE
    + fAliasSampler->SampleX<Backend>(deltaE,probNA,aliasInd,icol,fraction);
  sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
IonisationMoller::SampleSinTheta(typename Backend::Double_v energyIn,
                                 typename Backend::Double_v energyOut) const
{
  typedef Mask_v<typename Backend::Double_v>   Bool_t;
  using Double_v = typename Backend::Double_v;

  //angle of the scatterred electron

  Double_v energy = energyIn + electron_mass_c2;
  Double_v totalMomentum = Sqrt(energyIn*(energyIn + 2.0*electron_mass_c2));

  Double_v deltaMomentum = Sqrt(energyOut * (energyOut + 2.0*electron_mass_c2));
  Double_v cost =  energyOut * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum);

  Double_v sint2 = (1.0 - cost)*(1. + cost);
  Double_v sinTheta;
  Bool_t condition2 = sint2 < 0.0;
  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );

  return sinTheta;
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
IonisationMoller::InteractKernelCR(typename Backend::Double_v  kineticEnergy,
                                   Index_v<typename Backend::Double_v>   zElement,
                                   typename Backend::Double_v& deltaKinEnergy,
                                   typename Backend::Double_v& sinTheta)
{
  typedef Mask_v<typename Backend::Double_v> Bool_t;
  using Double_v = typename Backend::Double_v;

  //temporary - set by material
  Double_v cutEnergy = fDeltaRayThreshold;
  Double_v maxEnergy = 1.0*TeV;

  //based on G4MollerBhabhaModel::SampleSecondaries
  Double_v tmin = cutEnergy;
  Double_v tmax = 0.5*kineticEnergy;

  Bool_t condCut = (tmax < maxEnergy);
  MaskedAssign(!condCut, maxEnergy, &tmax);

  condCut |= (tmax >= tmin );

  if(Backend::early_returns && IsEmpty(condCut)) return;

  Double_v energy = kineticEnergy + electron_mass_c2;
  Double_v xmin   = tmin/kineticEnergy;
  Double_v xmax   = tmax/kineticEnergy;
  Double_v gam    = energy/electron_mass_c2;
  Double_v gamma2 = gam*gam;

  //Moller (e-e-) scattering
  Double_v gg = (2.0*gam - 1.0)/gamma2;

  Double_v x = SampleSequential<Backend>(xmin,xmax,gg);

  deltaKinEnergy = x * kineticEnergy;

  Double_v totalMomentum = Sqrt(kineticEnergy*(kineticEnergy + 2.0*electron_mass_c2));

  Double_v deltaMomentum =
    Sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  Double_v cost = deltaKinEnergy * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum );

  Bool_t condCos = (cost <= 1.0);
  MaskedAssign(!condCos, 1.0, &cost);

  Double_v sint2 = (1.0 - cost)*(1.0 + cost);

  Bool_t condSin2 = (sint2 >= 0.0);
  Double_v zero(0.0);
  CondAssign(condSin2, Sqrt(sint2), zero, &sinTheta);
}

template<class Backend>
inline
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
IonisationMoller::SampleSequential(typename Backend::Double_v xmin,
                                   typename Backend::Double_v xmax,
                                   typename Backend::Double_v gg) const
{
  typedef typename Backend::Int_t Int_t;
  using Double_v = typename Backend::Double_v;

  Double_v  q;
  Double_v  x;
  Double_v  z;

  Double_v  y = 1.0 - xmax;
  Double_v grej =  1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));
  do {
    q = UniformRandom<Backend>(fRandomState,Int_t(fThreadId));
    x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
    y = 1.0 - x;
    z = 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y));
  } while(grej * UniformRandom<Backend>(fRandomState,Int_t(fThreadId)) > z);

  return x;
}

#ifndef VECCORE_NVCC
template<>
inline
VECCORE_CUDA_HOST_DEVICE
typename kVc::Double_v
IonisationMoller::SampleSequential<kVc>(typename kVc::Double_v xmin,
                                        typename kVc::Double_v xmax,
                                        typename kVc::Double_v gg) const
{
  typedef typename kVc::Double_v Double_v;

  Double_v  x;
  Double_v  y = 1.0 - xmax;
  Double_v  grej =  1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));

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
VECCORE_CUDA_HOST_DEVICE void
IonisationMoller::InteractKernelUnpack(typename Backend::Double_v energyIn,
                                       Index_v<typename Backend::Double_v>   zElement,
                                       typename Backend::Double_v& energyOut,
                                       typename Backend::Double_v& sinTheta,
                                       Mask_v<typename Backend::Double_v> &status)
{
  //dummy for now
  energyOut = energyIn;
  sinTheta =  0;
}

} // end namespace impl
} // end namespace vecphys

#endif
