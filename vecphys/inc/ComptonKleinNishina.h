#ifndef ComptonKleinNishina_H
#define ComptonKleinNishina_H 1

#include "base/Global.h"
#include "base/PhysicalConstants.h"

#include "GUTrack.h"
#include "GUAliasSampler.h"

#include "EmModelBase.h"


namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class ComptonKleinNishina : public EmModelBase<ComptonKleinNishina>
{
public:

  VECCORE_CUDA_HOST
  ComptonKleinNishina(Random_t* states = 0, int threadId = -1);

  VECCORE_CUDA_HOST_DEVICE
  ComptonKleinNishina(Random_t* states, int threadId, GUAliasSampler* sampler);

  VECCORE_CUDA_HOST_DEVICE
    ~ComptonKleinNishina();//{}

  VECCORE_CUDA_HOST
  void Initialization();

  //interfaces for tables
  VECCORE_CUDA_HOST
  void BuildCrossSectionTablePerAtom(int Z);

  VECCORE_CUDA_HOST
  void BuildPdfTable(int Z, double *p);

public:
  // Auxiliary methods
  VECCORE_CUDA_HOST_DEVICE
  GUAliasSampler* GetSampler() {return fAliasSampler;}

  VECCORE_CUDA_HOST_DEVICE
  void SetSampler(GUAliasSampler* sampler) { fAliasSampler = sampler ;}

  //Alternative Interact method to test energy dependent subtasks for a
  //specific model. Eventually this method should replace the Interact
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
  VECCORE_CUDA_HOST_DEVICE typename
  Backend::Double_v
  CrossSectionKernel(typename Backend::Double_v energyIn,
                     Index_v<typename Backend::Double_v> zElement);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE void
  InteractKernel(typename Backend::Double_v energyIn,
                 Index_v<typename Backend::Double_v> zElement,
                 typename Backend::Double_v& energyOut,
                 typename Backend::Double_v& sinTheta);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE void
  InteractKernelCR(typename Backend::Double_v energyIn,
                   Index_v<typename Backend::Double_v> zElement,
                   typename Backend::Double_v& energyOut,
                   typename Backend::Double_v& sinTheta);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE void
  InteractKernelUnpack(typename Backend::Double_v energyIn,
                       Index_v<typename Backend::Double_v> zElement,
                       typename Backend::Double_v& energyOut,
                       typename Backend::Double_v& sinTheta,
                       Mask_v<typename Backend::Double_v> &status);

  template<class Backend>
  inline
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  SampleSequential(typename Backend::Double_v E0_m,
                   typename Backend::Double_v test,
                   typename Backend::Double_v alpha1,
                   typename Backend::Double_v epsil0sq,
                   typename Backend::Double_v &sint2);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  SampleSinTheta(typename Backend::Double_v energyIn,
                 typename Backend::Double_v energyOut);

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

  VECCORE_CUDA_HOST_DEVICE
  GUTrack& GetSecondaryElectron() { return fSecondaryElectron; }

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<ComptonKleinNishina>;

private:

  GUTrack fSecondaryElectron;
};

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
ComptonKleinNishina::CrossSectionKernel(typename Backend::Double_v energy,
                                        Index_v<typename Backend::Double_v> Z)
{
  using Double_v = typename Backend::Double_v;

  Double_v Z2 = Z*Z;
  Double_v p1 =  2.7965e-1 +  1.9756e-5*Z + -3.9178e-7*Z2;
  Double_v p2 = -1.8300e-1 + -1.0205e-2*Z +  6.8241e-5*Z2;
  Double_v p3 =  6.7527    + -7.3913e-2*Z +  6.0480e-5*Z2;
  Double_v p4 = -1.9798e+1 +  2.7079e-2*Z +  3.0274e-4*Z2;

  Mask_v<Double_v> condZ = (Z < 1.5);
  Double_v T0 = Blend(condZ, Double_v(15.*keV), Double_v(40.*keV));

  Double_v X  =  math::Max(energy,T0)/electron_mass_c2;
  Double_v X2 = X*X;
  Double_v sigma = p1*math::Log(1.+2.*X)/X
          + (p2 + p3*X + p4*X2)/(1. + 20.*X + 230.*X2 + 440.*X2*X);

  Double_v sigmaOut = Z*sigma*barn;

  Mask_v<Double_v> condE = (energy > T0);

  if(EarlyReturnAllowed() && MaskFull(condE)) return sigmaOut;

  //correction when energy < T0
  Double_v dT0 = 1.*keV;
  X = (T0+dT0) / electron_mass_c2 ;
  sigma = p1*math::Log(1.+2.*X)/X
          + (p2 + p3*X + p4*X2)/(1. + 20.*X + 230.*X2 + 440.*X2*X);

  Double_v c1 = -T0*(Z*sigma*barn-sigmaOut)/(sigmaOut*dT0);
  Double_v c2 = Blend(condZ, Double_v(0.150),
                             Double_v(0.375) - Double_v(0.0556) * math::Log(Double_v(Z)));
  Double_v  y = math::Log(energy/T0);

  MaskedAssign(sigmaOut, !condE, sigmaOut*math::Exp(-y*(c1+c2*y)));
  MaskedAssign(sigmaOut, energy < fLowEnergyLimit, 0.0);

  return  sigmaOut;
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
ComptonKleinNishina::InteractKernel(typename Backend::Double_v  energyIn,
                                    Index_v<typename Backend::Double_v> zElement,
                                    typename Backend::Double_v& energyOut,
                                    typename Backend::Double_v& sinTheta)
{
  using Double_v = typename Backend::Double_v;

  Index_v<Double_v> irow;
  Index_v<Double_v> icol;
  Double_v fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,irow,icol,fraction);

  Double_v probNA;
  Index_v<Double_v> aliasInd;

  Double_v ncol(fAliasSampler->GetSamplesPerEntry());
  Index_v<Double_v> index = ncol*irow + icol;
  fAliasSampler->GatherAlias<Backend>(index,probNA,aliasInd);

  Double_v mininumE = energyIn/(1+2.0*energyIn*inv_electron_mass_c2);
  Double_v deltaE = energyIn - mininumE;

  energyOut = mininumE + fAliasSampler->SampleXL<Backend>(zElement,
                                        deltaE,probNA,aliasInd,irow,icol);
  sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
ComptonKleinNishina::InteractKernelCR(typename Backend::Double_v energyIn,
                                      Index_v<typename Backend::Double_v> zElement,
                                      typename Backend::Double_v& energyOut,
                                      typename Backend::Double_v& sinTheta)
{
  using Double_v = typename Backend::Double_v;

  Double_v E0_m = energyIn/electron_mass_c2;

  Double_v eps0 = 1./(1. + 2.*E0_m);
  Double_v epsilon0sq = eps0*eps0;
  Double_v alpha1     = - math::Log(eps0);
  Double_v alpha2  = 0.5*(1.- epsilon0sq);

  Double_v test = alpha1/(alpha1+alpha2);

  Double_v sint2;
  Double_v epsilon = SampleSequential<Backend>(E0_m,test,alpha1,epsilon0sq,sint2);

  energyOut = epsilon*energyIn;
  sinTheta = math::Sqrt(sint2);
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
ComptonKleinNishina::SampleSinTheta(typename Backend::Double_v energyIn,
                                    typename Backend::Double_v energyOut)
{
  using Double_v = typename Backend::Double_v;

  //angle of the scatterred photon

  Double_v epsilon = Blend(energyOut > energyIn, Double_v(1.0), energyOut/energyIn);
  Double_v E0_m    = inv_electron_mass_c2*energyIn;
  Double_v onecost = (1.0 - epsilon)/(epsilon*E0_m);
  Double_v sint2   = onecost*(2.-onecost);

  return Blend(sint2 < 0.0, Double_v(0.0), math::Sqrt(sint2));
}

template<class Backend>
inline
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
ComptonKleinNishina::SampleSequential(typename Backend::Double_v E0_m,
                                      typename Backend::Double_v test,
                                      typename Backend::Double_v alpha1,
                                      typename Backend::Double_v epsil0sq,
                                      typename Backend::Double_v &sint2)
{
  using Double_v = typename Backend::Double_v;

  Double_v epsilon;
  Double_v greject;
  Mask_v<Double_v> done(false);

  do {
    Double_v tmp = Blend(test > UniformRandom<Double_v>(&fRandomState, &fThreadId),
                         math::Exp(-alpha1*UniformRandom<Double_v>(&fRandomState, &fThreadId)),
                         math::Sqrt(epsil0sq+(1.0 - epsil0sq)*UniformRandom<Double_v>(&fRandomState, &fThreadId)));
    MaskedAssign(epsilon, !done, tmp);
    Double_v onecost = (1.0 - epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);
    done |= greject > UniformRandom<Double_v>(&fRandomState, &fThreadId);
  } while (!MaskFull(done));

  return epsilon;
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
ComptonKleinNishina::InteractKernelUnpack(typename Backend::Double_v energyIn,
                                          Index_v<typename Backend::Double_v> zElement,
                                          typename Backend::Double_v& energyOut,
                                          typename Backend::Double_v& sinTheta,
                                          Mask_v<typename Backend::Double_v>& status)
{
  using Double_v = typename Backend::Double_v;

  Double_v E0_m = energyIn/electron_mass_c2;

  Double_v eps0 = 1./(1. + 2.*E0_m);
  Double_v epsilon0sq = eps0*eps0;
  Double_v alpha1     = - math::Log(eps0);
  Double_v alpha2  = 0.5*(1.- epsilon0sq);

  Double_v test = alpha1/(alpha1+alpha2);

  Double_v epsilon = Blend(test > UniformRandom<Double_v>(&fRandomState, &fThreadId),
                           math::Exp(-alpha1*UniformRandom<Double_v>(&fRandomState, &fThreadId)),
                           math::Sqrt(epsilon0sq+(1.- epsilon0sq)*UniformRandom<Double_v>(&fRandomState, &fThreadId)));

  Double_v onecost = (1.- epsilon)/(epsilon*E0_m);
  Double_v sint2   = onecost*(2.-onecost);

  Double_v greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);

  status = greject < UniformRandom<Double_v>(&fRandomState, &fThreadId);

  energyOut = epsilon*energyIn;
  sinTheta = math::Sqrt(sint2);
}

//Alternative Interact method

template <typename Backend>
VECCORE_CUDA_HOST_DEVICE
void ComptonKleinNishina::ModelInteract(GUTrack&  inProjectile,
                                        const int targetElement,
                                        GUTrack&  outSecondary )
{
  double energyIn = inProjectile.E;

  //check for the validity of energy
  if(energyIn < fLowEnergyLimit || energyIn > fHighEnergyLimit) return;

  double energyOut =0;
  double sinTheta = 0;

  const double aliaslimit = 100.0*MeV;

  if(energyIn< aliaslimit) {
    InteractKernel<Backend>(energyIn,targetElement,energyOut,sinTheta);
  }
  else {
    InteractKernelCR<Backend>(energyIn,targetElement,energyOut,sinTheta);
  }

  //update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn,energyOut,sinTheta,
                                inProjectile,outSecondary);
}

#ifndef VECCORE_NVCC
template <typename Backend>
void ComptonKleinNishina::ModelInteract(GUTrack_v& inProjectile,
                                        const int* targetElements,
                                        GUTrack_v& outSecondary)
{
  using Double_v = typename Backend::Double_v;

  //check for the validity of energy
  int nTracks = inProjectile.numTracks;

  // this inclusive check may be redundant as this model/process should not be
  // selected if energy of the track is outside the valid energy region
  //  if(inProjectile.E[0]         < fLowEnergyLimit ||
  //     inProjectile.E[nTracks-1] > fHighEnergyLimit) return;

  //filtering the energy region for the alias method - setable if necessary
  const double aliaslimit = 100.0*MeV;
  double* start = inProjectile.E;
  auto indexAliasLimit = std::lower_bound(start,start+nTracks,aliaslimit) - start;

  for(int j = 0; j < nTracks  ; ++j) {
    assert( (targetElements[j] > 0)  && (targetElements[j] <= maximumZ ) );
  }

  int ibase= 0;
  int numChunks= (nTracks/VectorSize<Double_v>());

  for(int i= 0; i < numChunks ; ++i) {

    Double_v energyIn(&inProjectile.E[ibase]);
    Double_v sinTheta(0.);
    Double_v energyOut;

    Index_v<Double_v> zElement(targetElements[ibase]);

    if(ibase < indexAliasLimit) {
      InteractKernel<Backend>(energyIn,zElement,energyOut,sinTheta);
    }
    else {
      InteractKernelCR<Backend>(energyIn,zElement,energyOut,sinTheta);
    }

    ConvertXtoFinalState<Backend>(energyIn, energyOut, sinTheta, ibase, inProjectile, outSecondary);

    ibase+= VectorSize<Double_v>();
  }

  //leftover - do scalar (temporary)
  for(int i = numChunks*VectorSize<Double_v>(); i < inProjectile.numTracks; ++i) {

    double senergyIn= inProjectile.E[i];
    double senergyOut, ssinTheta;
    //use InteractKernel for any leftover to be consistent with EmBaseModel
    InteractKernel<backend::Scalar>(senergyIn,targetElements[i],senergyOut,ssinTheta);
    ConvertXtoFinalState_Scalar<backend::Scalar>(senergyIn, senergyOut, ssinTheta, i, inProjectile, outSecondary);
  }
}

#endif

} // end namespace impl
} // end namespace vecphys

#endif