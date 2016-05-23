#ifndef IonisationMoller_H
#define IonisationMoller_H 1

#include "base/VPGlobal.h"
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
                 typename Backend::Double_v energyOut);

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
                   typename Backend::Double_v gg);

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

  using Double_v = typename Backend::Double_v;

  Double_v sigmaOut = 0.;
  Mask_v<Double_v> belowLimit = Mask_v<Double_v>(false);
  //low energy limit
  belowLimit |= ( energy < fLowEnergyLimit );
  if(EarlyReturnAllowed() && MaskFull(belowLimit)) return;

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
          - gg*math::Log( xmax*(1.0 - xmin)/(xmin*(1.0 - xmax)) ) ) / beta2;

  sigmaOut *= Z*twopi_mc2_rcl2/energy;

  //this is the case if one of E < belowLimit
  MaskedAssign(sigmaOut, belowLimit, Double_v(0.0));
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
IonisationMoller::InteractKernel(typename Backend::Double_v  energyIn,
                                 Index_v<typename Backend::Double_v>   zElement,
                                 typename Backend::Double_v& energyOut,
                                 typename Backend::Double_v& sinTheta)
{
  using Double_v = typename Backend::Double_v;

  Index_v<Double_v>   irow;
  Index_v<Double_v>   icol;
  Double_v  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,irow,icol,fraction);

  Double_v probNA;
  Index_v<Double_v>  aliasInd;

  //this did not used to work - Fixed SW
  Double_v ncol(fAliasSampler->GetSamplesPerEntry());
  Index_v<Double_v>   index = ncol*irow + icol;
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
                                 typename Backend::Double_v energyOut)
{
  using Double_v = typename Backend::Double_v;

  //angle of the scatterred electron

  Double_v energy = energyIn + electron_mass_c2;
  Double_v totalMomentum = math::Sqrt(energyIn*(energyIn + 2.0*electron_mass_c2));
  Double_v deltaMomentum = math::Sqrt(energyOut * (energyOut + 2.0*electron_mass_c2));
  Double_v cost =  energyOut * (energy + electron_mass_c2) / (deltaMomentum * totalMomentum);
  Double_v sint2 = 1.0 - cost*cost;

  return Blend(sint2 < 0.0, Double_v(0.0), math::Sqrt(sint2));
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
IonisationMoller::InteractKernelCR(typename Backend::Double_v  kineticEnergy,
                                   Index_v<typename Backend::Double_v>   zElement,
                                   typename Backend::Double_v& deltaKinEnergy,
                                   typename Backend::Double_v& sinTheta)
{
  using Double_v = typename Backend::Double_v;

  //temporary - set by material
  Double_v cutEnergy = fDeltaRayThreshold;
  Double_v maxEnergy = 1.0*TeV;

  //based on G4MollerBhabhaModel::SampleSecondaries
  Double_v tmin = cutEnergy;
  Double_v tmax = 0.5*kineticEnergy;

  Mask_v<Double_v> condCut = (tmax < maxEnergy);
  MaskedAssign(tmax, !condCut, maxEnergy);

  condCut |= (tmax >= tmin );

  if(EarlyReturnAllowed() && MaskEmpty(condCut)) return;

  Double_v energy = kineticEnergy + electron_mass_c2;
  Double_v xmin   = tmin/kineticEnergy;
  Double_v xmax   = tmax/kineticEnergy;
  Double_v gam    = energy/electron_mass_c2;
  Double_v gamma2 = gam*gam;

  //Moller (e-e-) scattering
  Double_v gg = (2.0*gam - 1.0)/gamma2;

  Double_v x = SampleSequential<Backend>(xmin,xmax,gg);

  deltaKinEnergy = x * kineticEnergy;

  Double_v totalMomentum = math::Sqrt(kineticEnergy*(kineticEnergy + 2.0*electron_mass_c2));

  Double_v deltaMomentum =
    math::Sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  Double_v cost = deltaKinEnergy * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum );

  Mask_v<Double_v> condCos = (cost <= 1.0);
  MaskedAssign(cost, !condCos, Double_v(1.0));

  Double_v sint2 = (1.0 - cost)*(1.0 + cost);

  Mask_v<Double_v> condSin2 = (sint2 >= 0.0);
  sinTheta = Blend(condSin2, math::Sqrt(sint2), Double_v(0.0));
}

template<class Backend>
inline
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
IonisationMoller::SampleSequential(typename Backend::Double_v xmin,
                                   typename Backend::Double_v xmax,
                                   typename Backend::Double_v gg)
{
  using Double_v = typename Backend::Double_v;

  Double_v  q;
  Double_v  x;
  Double_v  z;
  Mask_v<Double_v> done(false);

  Double_v  y = 1.0 - xmax;
  Double_v grej =  1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));
  do {
    q = UniformRandom<Double_v>(&fRandomState, &fThreadId);
    MaskedAssign(x, !done, xmin*xmax/(xmin*(1.0 - q) + xmax*q));
    MaskedAssign(y, !done, 1.0 - x);
    MaskedAssign(z, !done, 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y)));
    done |= z < grej * UniformRandom<Double_v>(&fRandomState, &fThreadId);
  } while(!MaskFull(done));

  return x;
}

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
