#ifndef ConversionBetheHeitler_H
#define ConversionBetheHeitler_H 1

#include "backend/Backend.h"
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

  VECPHYS_CUDA_HEADER_HOST
  ConversionBetheHeitler(Random_t* states = 0, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  ConversionBetheHeitler(Random_t* states, int threadId, 
                         GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH
  ~ConversionBetheHeitler(){}

  VECPHYS_CUDA_HEADER_HOST
  void Initialization();

  //interfaces for tables
  VECPHYS_CUDA_HEADER_HOST 
  void BuildCrossSectionTablePerAtom(int Z);

  VECPHYS_CUDA_HEADER_HOST
  void BuildPdfTable(int Z, double *p);

public:
  // Implementation methods
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH 
  typename Backend::Double_t
  CrossSectionKernel(typename Backend::Double_t  energyIn,
                     typename Backend::Index_t   zElement);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void 
  InteractKernel(typename Backend::Double_t  energyIn, 
                 typename Backend::Index_t   zElement,
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
  void
  SampleSinTheta(typename Backend::Double_t energyElectron,
                 typename Backend::Double_t energyPositron,
		 typename Backend::Double_t& sinThetaElectron,
		 typename Backend::Double_t& sinThetaPositron) const; 

  VECPHYS_CUDA_HEADER_BOTH 
  void SampleByCompositionRejection(int    elementZ,
                                    double energyIn,
                                    double& energyOut,
                                    double& sinTheta);

  VECPHYS_CUDA_HEADER_BOTH 
  double 
  GetG4CrossSection(double  energyIn, 
                    const int zElement);

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection(int Zelement, 
                                   double Ein,
                                   double outEphoton );

  //this should be a method of GUElement 
  /*
  VECPHYS_CUDA_HEADER_BOTH 
  double ComputeCoulombFactor(double Zeff) const;
  */
  VECPHYS_CUDA_HEADER_BOTH 
  double ScreenFunction1(double screenVariable) const;

  VECPHYS_CUDA_HEADER_BOTH 
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
VECPHYS_CUDA_HEADER_BOTH 
void 
ConversionBetheHeitler::InteractKernel(typename Backend::Double_t  energyIn, 
                                       typename Backend::Index_t   zElement,
                                       typename Backend::Double_t& energyOut,
                                       typename Backend::Double_t& sinTheta)
{
  // now return only secondary electron information and
  // a positron will be created based on the electron - eventually we need a common
  // interface  to fill produced secondaries into a single stact
  typedef typename Backend::Bool_t   Bool_t;
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
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);
  
  Double_t mininumE = electron_mass_c2;
  Double_t deltaE = energyIn - mininumE;

  //electron energy
  energyOut = mininumE + fAliasSampler->SampleX<Backend>(deltaE,probNA,
					        aliasInd,icol,fraction);

  Double_t r1 = UniformRandom<Backend>(fRandomState,fThreadId);
  Bool_t condition = 0.5 > r1;

  Double_t energyElectron;
  Double_t energyPositron;

  //check correctness
  MaskedAssign( condition, energyOut, &energyElectron);     
  MaskedAssign( condition, energyIn - energyOut, &energyPositron);     

  MaskedAssign(!condition, energyOut, &energyPositron);     
  MaskedAssign(!condition, energyIn - energyOut, &energyElectron);     

  Double_t sinThetaElectron;
  Double_t sinThetaPositron;
  SampleSinTheta<Backend>(energyElectron,energyPositron,
                          sinThetaElectron, sinThetaPositron);

  //fill secondaries
  energyOut = energyElectron;
  sinTheta = sinThetaElectron;
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void
ConversionBetheHeitler::
SampleSinTheta(typename Backend::Double_t energyElectron,
               typename Backend::Double_t energyPositron,
	       typename Backend::Double_t& sinThetaElectron,
	       typename Backend::Double_t& sinThetaPositron) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  //angles of the pair production (gamma -> e+e-)

  Double_t u;
  const double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  Double_t r1 =  UniformRandom<Backend>(fRandomState,fThreadId);
  Bool_t condition = 9./(9. + d) > r1;
  MaskedAssign( condition, -log( UniformRandom<Backend>(fRandomState,fThreadId)*
                       UniformRandom<Backend>(fRandomState,fThreadId))/a1, &u );
  MaskedAssign(!condition, -log( UniformRandom<Backend>(fRandomState,fThreadId)*
                       UniformRandom<Backend>(fRandomState,fThreadId))/a2, &u );

  Double_t TetEl = u*electron_mass_c2/energyElectron;
  Double_t TetPo = u*electron_mass_c2/energyPositron;

  //sinTheta - just return theta instead!
  sinThetaElectron =  sin(TetEl);
  sinThetaPositron = -sin(TetPo);
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t 
ConversionBetheHeitler::
CrossSectionKernel(typename Backend::Double_t energy,
                   typename Backend::Index_t Z)
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  Double_t sigma = 0.;

  if ( Z < 0.9 || energy <= 2.0*electron_mass_c2 ) { return sigma; }

  Double_t energySave = energy;

  //gamma energyLimit = 1.5*MeV
  Double_t energyLimit = 1.5*MeV;
  Bool_t condition = energy < energyLimit;
  MaskedAssign( condition, energyLimit, &energy );
  
  Double_t X = log(energy/electron_mass_c2);
  Double_t X2 = X*X;
  Double_t X3 =X2*X;
  Double_t X4 =X3*X;
  Double_t X5 =X4*X;

  //put coff's to a constant header
  /*
  Double_t a0= 8.7842e+2*microbarn; 
  Double_t a1=-1.9625e+3*microbarn; 
  Double_t a2= 1.2949e+3*microbarn;
  Double_t a3=-2.0028e+2*microbarn; 
  Double_t a4= 1.2575e+1*microbarn; 
  Double_t a5=-2.8333e-1*microbarn;
  
  Double_t b0=-1.0342e+1*microbarn; 
  Double_t b1= 1.7692e+1*microbarn; 
  Double_t b2=-8.2381   *microbarn;
  Double_t b3= 1.3063   *microbarn; 
  Double_t b4=-9.0815e-2*microbarn; 
  Double_t b5= 2.3586e-3*microbarn;
  
  Double_t c0=-4.5263e+2*microbarn; 
  Double_t c1= 1.1161e+3*microbarn; 
  Double_t c2=-8.6749e+2*microbarn;
  Double_t c3= 2.1773e+2*microbarn; 
  Double_t c4=-2.0467e+1*microbarn; 
  Double_t c5= 6.5372e-1*microbarn;
  */

  Double_t F1 = a0 + a1*X + a2*X2 + a3*X3 + a4*X4 + a5*X5;
  Double_t F2 = b0 + b1*X + b2*X2 + b3*X3 + b4*X4 + b5*X5;
  Double_t F3 = c0 + c1*X + c2*X2 + c3*X3 + c4*X4 + c5*X5;     

  sigma = (Z + 1.)*(F1*Z + F2*Z*Z + F3);
  Bool_t done = energySave < energyLimit;

  if(Any(done)) {
    X = (energySave - 2.*electron_mass_c2)/(energyLimit - 2.*electron_mass_c2);
    Double_t tmpsigma = sigma*X*X;
    MaskedAssign( done, tmpsigma, &sigma );
  }

  Bool_t check = sigma < 0.;
  MaskedAssign( check, 0., &sigma );

  return sigma*microbarn;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
ConversionBetheHeitler::InteractKernelCR(typename Backend::Double_t  energyIn, 
                                         typename Backend::Index_t   zElement,
                                         typename Backend::Double_t& energyOut,
                                         typename Backend::Double_t& sinTheta)
{
  //dummy for now
  energyOut = 0.0;
  sinTheta = 0.0;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
ConversionBetheHeitler::InteractKernelUnpack(typename Backend::Double_t energyIn, 
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
