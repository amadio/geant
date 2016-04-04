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
                       typename Backend::Bool_t &status);

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  void
  SampleSinTheta(typename Backend::Double_v energyElectron,
                 typename Backend::Double_v energyPositron,
		 typename Backend::Double_v& sinThetaElectron,
		 typename Backend::Double_v& sinThetaPositron) const;

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
  typedef typename Backend::Bool_t   Bool_t;
  typedef Index_v<typename Backend::Double_v>  Index_t;
  using Double_v = typename Backend::Double_v;

  //early return if E_gamma < 2*electron_mass_c2

  Index_t   irow;
  Index_t   icol;
  Double_v  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,irow,icol,fraction);

  Double_v probNA;
  Index_t  aliasInd;

  //this did not used to work - Fixed SW
  Double_v ncol(fAliasSampler->GetSamplesPerEntry());
  Index_t   index = ncol*irow + icol;
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);

  Double_v mininumE = electron_mass_c2;
  Double_v deltaE = energyIn - mininumE;

  //electron energy
  energyOut = mininumE + fAliasSampler->SampleX<Backend>(deltaE,probNA,
					        aliasInd,icol,fraction);

  Double_v r1 = UniformRandom<Backend>(fRandomState,fThreadId);
  Bool_t condition = 0.5 > r1;

  Double_v energyElectron;
  Double_v energyPositron;

  //check correctness
  MaskedAssign( condition, energyOut, &energyElectron);
  MaskedAssign( condition, energyIn - energyOut, &energyPositron);

  MaskedAssign(!condition, energyOut, &energyPositron);
  MaskedAssign(!condition, energyIn - energyOut, &energyElectron);

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
	       typename Backend::Double_v& sinThetaPositron) const
{
  typedef typename Backend::Bool_t   Bool_t;
  using Double_v = typename Backend::Double_v;

  //angles of the pair production (gamma -> e+e-)

  Double_v u;
  const double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  Double_v r1 =  UniformRandom<Backend>(fRandomState,fThreadId);
  Bool_t condition = 9./(9. + d) > r1;
  MaskedAssign( condition, -log( UniformRandom<Backend>(fRandomState,fThreadId)*
                       UniformRandom<Backend>(fRandomState,fThreadId))/a1, &u );
  MaskedAssign(!condition, -log( UniformRandom<Backend>(fRandomState,fThreadId)*
                       UniformRandom<Backend>(fRandomState,fThreadId))/a2, &u );

  Double_v TetEl = u*electron_mass_c2/energyElectron;
  Double_v TetPo = u*electron_mass_c2/energyPositron;

  //sinTheta - just return theta instead!
  sinThetaElectron =  sin(TetEl);
  sinThetaPositron = -sin(TetPo);
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
ConversionBetheHeitler::
CrossSectionKernel(typename Backend::Double_v energy,
                   Index_v<typename Backend::Double_v> Z)
{
  typedef typename Backend::Bool_t   Bool_t;
  using Double_v = typename Backend::Double_v;

  Double_v sigma = 0.;

  if ( Z < 0.9 || energy <= 2.0*electron_mass_c2 ) { return sigma; }

  Double_v energySave = energy;

  //gamma energyLimit = 1.5*MeV
  Double_v energyLimit = 1.5*MeV;
  Bool_t condition = energy < energyLimit;
  MaskedAssign( condition, energyLimit, &energy );

  Double_v X = log(energy/electron_mass_c2);
  Double_v X2 = X*X;
  Double_v X3 =X2*X;
  Double_v X4 =X3*X;
  Double_v X5 =X4*X;

  //put coff's to a constant header
  /*
  Double_v a0= 8.7842e+2*microbarn;
  Double_v a1=-1.9625e+3*microbarn;
  Double_v a2= 1.2949e+3*microbarn;
  Double_v a3=-2.0028e+2*microbarn;
  Double_v a4= 1.2575e+1*microbarn;
  Double_v a5=-2.8333e-1*microbarn;

  Double_v b0=-1.0342e+1*microbarn;
  Double_v b1= 1.7692e+1*microbarn;
  Double_v b2=-8.2381   *microbarn;
  Double_v b3= 1.3063   *microbarn;
  Double_v b4=-9.0815e-2*microbarn;
  Double_v b5= 2.3586e-3*microbarn;

  Double_v c0=-4.5263e+2*microbarn;
  Double_v c1= 1.1161e+3*microbarn;
  Double_v c2=-8.6749e+2*microbarn;
  Double_v c3= 2.1773e+2*microbarn;
  Double_v c4=-2.0467e+1*microbarn;
  Double_v c5= 6.5372e-1*microbarn;
  */

  Double_v F1 = a0 + a1*X + a2*X2 + a3*X3 + a4*X4 + a5*X5;
  Double_v F2 = b0 + b1*X + b2*X2 + b3*X3 + b4*X4 + b5*X5;
  Double_v F3 = c0 + c1*X + c2*X2 + c3*X3 + c4*X4 + c5*X5;

  sigma = (Z + 1.)*(F1*Z + F2*Z*Z + F3);
  Bool_t done = energySave < energyLimit;

  if(Any(done)) {
    X = (energySave - 2.*electron_mass_c2)/(energyLimit - 2.*electron_mass_c2);
    Double_v tmpsigma = sigma*X*X;
    MaskedAssign( done, tmpsigma, &sigma );
  }

  Bool_t check = sigma < 0.;
  MaskedAssign( check, 0., &sigma );

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
                                             typename Backend::Bool_t &status)
{
  //dummy for now
  energyOut = energyIn;
  sinTheta =  0;
}

} // end namespace impl
} // end namespace vecphys

#endif
