#ifndef ConversionBetheHeitler_H
#define ConversionBetheHeitler_H 1

#include "base/Global.h"
#include "base/PhysicalConstants.h"

// #include "GUAuxFunctions.h"    // Define sincos if needed

#include "GUConstants.h"
#include "GUTrack.h"

#include "EmModelBase.h"


namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class ConversionBetheHeitler : public EmModelBase<ConversionBetheHeitler>
{
public:

  VECCORE_CUDA_HOST
  ConversionBetheHeitler(Random_t* states = 0, int threadId = -1);

  VECCORE_CUDA_HOST_DEVICE
  ConversionBetheHeitler(Random_t* states, int threadId,
                         GUAliasSampler* sampler);

  VECCORE_CUDA_HOST_DEVICE
  ~ConversionBetheHeitler(){}

  VECCORE_CUDA_HOST
  void Initialization();

  //interfaces for tables
  VECCORE_CUDA_HOST
  void BuildCrossSectionTablePerAtom(int Z);

  VECCORE_CUDA_HOST
  void BuildPdfTable(int Z, double *p);

public:
  // Implementation methods
  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  CrossSectionKernel(typename Backend::Double_v  energyIn,
                     Index_v<typename Backend::Double_v>   zElement);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE void
  InteractKernel(typename Backend::Double_v  energyIn,
                 Index_v<typename Backend::Double_v>   zElement,
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
  void
  SampleSinTheta(typename Backend::Double_v energyElectron,
                 typename Backend::Double_v energyPositron,
		 typename Backend::Double_v& sinThetaElectron,
		 typename Backend::Double_v& sinThetaPositron);

  VECCORE_CUDA_HOST_DEVICE
  void SampleByCompositionRejection(int    elementZ,
                                    double energyIn,
                                    double& energyOut,
                                    double& sinTheta);

  VECCORE_CUDA_HOST_DEVICE
  double
  GetG4CrossSection(double  energyIn,
                    const int zElement);

  VECCORE_CUDA_HOST_DEVICE
  double CalculateDiffCrossSection(int Zelement,
                                   double Ein,
                                   double outEphoton );

  //this should be a method of GUElement
  /*
  VECCORE_CUDA_HOST_DEVICE
  double ComputeCoulombFactor(double Zeff) const;
  */
  VECCORE_CUDA_HOST_DEVICE
  double ScreenFunction1(double screenVariable) const;

  VECCORE_CUDA_HOST_DEVICE
  double ScreenFunction2(double screenVariable) const;

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<ConversionBetheHeitler>;

private:
  //Secondary
  GUTrack fElectron;
  GUTrack fPositron;
};

//Implementation

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
void
ConversionBetheHeitler::InteractKernel(typename Backend::Double_v  energyIn,
                                       Index_v<typename Backend::Double_v>   zElement,
                                       typename Backend::Double_v& energyOut,
                                       typename Backend::Double_v& sinTheta)
{
  // now return only secondary electron information and
  // a positron will be created based on the electron - eventually we need a common
  // interface  to fill produced secondaries into a single stact
  using Double_v = typename Backend::Double_v;

  //early return if E_gamma < 2*electron_mass_c2

  Index_v<Double_v>   irow;
  Index_v<Double_v>   icol;
  Double_v  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,irow,icol,fraction);

  Double_v probNA;
  Index_v<Double_v>  aliasInd;

  //this did not used to work - Fixed SW
  Double_v ncol(fAliasSampler->GetSamplesPerEntry());
  Index_v<Double_v>   index = ncol*irow + icol;
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);

  Double_v mininumE = electron_mass_c2;
  Double_v deltaE = energyIn - mininumE;

  //electron energy
  energyOut = mininumE + fAliasSampler->SampleX<Backend>(deltaE,probNA,
					        aliasInd,icol,fraction);

  Double_v r1 = UniformRandom<Double_v>(&fRandomState, &fThreadId);
  Mask_v<Double_v> condition = 0.5 > r1;

  Double_v energyElectron = Blend(condition, energyOut, energyIn - energyOut);
  Double_v energyPositron = Blend(condition, energyIn - energyOut, energyOut);

  Double_v sinThetaElectron;
  Double_v sinThetaPositron;
  SampleSinTheta<Backend>(energyElectron,energyPositron,
                          sinThetaElectron, sinThetaPositron);

  //fill secondaries
  energyOut = energyElectron;
  sinTheta = sinThetaElectron;
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
void
ConversionBetheHeitler::
SampleSinTheta(typename Backend::Double_v energyElectron,
               typename Backend::Double_v energyPositron,
               typename Backend::Double_v& sinThetaElectron,
               typename Backend::Double_v& sinThetaPositron)
{
  using Double_v = typename Backend::Double_v;

  // angles of the pair production (gamma -> e+e-)

  Double_v r1 = UniformRandom<Double_v>(&fRandomState, &fThreadId);

  Double_v r2 = UniformRandom<Double_v>(&fRandomState, &fThreadId) *
                UniformRandom<Double_v>(&fRandomState, &fThreadId);

  Double_v u = -math::Log(r2) * Double_v(electron_mass_c2) *
                Blend(r1 < 0.25, Double_v(1.6), Double_v(1.6 / 3.0));

  Double_v TetEl = u / energyElectron;
  Double_v TetPo = u / energyPositron;

  sinThetaElectron =  math::Sin(TetEl);
  sinThetaPositron = -math::Sin(TetPo);
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
ConversionBetheHeitler::
CrossSectionKernel(typename Backend::Double_v energy,
                   Index_v<typename Backend::Double_v> Z)
{
  using Double_v = typename Backend::Double_v;

  if (MaskFull(Z < 0.9 || energy <= Double_v(2.0*electron_mass_c2)))
    return 0.0;

  Double_v energyLimit = 1.5*MeV;
  Double_v energySave = energy;
  Mask_v<Double_v> done = energy < energyLimit;
  energy = math::Min(energy, energyLimit);

  Double_v X  = -math::Log(energy/electron_mass_c2);
  Double_v X2 =  X*X;
  Double_v X3 = X2*X;
  Double_v X4 = X3*X;
  Double_v X5 = X4*X;

  Double_v F1 = a0 + a1*X + a2*X2 + a3*X3 + a4*X4 + a5*X5;
  Double_v F2 = b0 + b1*X + b2*X2 + b3*X3 + b4*X4 + b5*X5;
  Double_v F3 = c0 + c1*X + c2*X2 + c3*X3 + c4*X4 + c5*X5;

  Double_v sigma = (Z + 1.0)*(F1*Z + F2*Z*Z + F3);

  if(!MaskEmpty(done)) {
    X = (energySave - 2.*electron_mass_c2)/(energyLimit - 2.*electron_mass_c2);
    MaskedAssign(sigma, done, sigma*X*X);
  }

  MaskedAssign(sigma, sigma < 0.0, 0.0);

  return sigma*microbarn;
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
ConversionBetheHeitler::InteractKernelCR(typename Backend::Double_v  energyIn,
                                         Index_v<typename Backend::Double_v>   zElement,
                                         typename Backend::Double_v& energyOut,
                                         typename Backend::Double_v& sinTheta)
{
  //dummy for now
  energyOut = 0.0;
  sinTheta = 0.0;
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE void
ConversionBetheHeitler::InteractKernelUnpack(typename Backend::Double_v energyIn,
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
