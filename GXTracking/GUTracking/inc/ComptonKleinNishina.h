#ifndef ComptonKleinNishina_H
#define ComptonKleinNishina_H 1

#include "backend/Backend.h"

#include "PhysicalConstants.h"
#include "GUTrack.h"

#include "GUAliasSampler.h"

#include "EmModelBase.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class ComptonKleinNishina : public EmModelBase<ComptonKleinNishina>
{
public:

  VECPHYS_CUDA_HEADER_HOST
  ComptonKleinNishina(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  ComptonKleinNishina(Random_t* states, int threadId, 
                        GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH
  ~ComptonKleinNishina();

  // Initializes this class and its sampler 
  VECPHYS_CUDA_HEADER_HOST
  void BuildOneTable( int Z,
                      const double xmin,
                      const double xmax,
                      const int nrow,
                      const int ncol);

  VECPHYS_CUDA_HEADER_HOST
  void BuildPdfTable(int Z,
                     const double xmin,
                     const double xmax,
                     const int nrow,
                     const int ncol,
                     double *p);

  VECPHYS_CUDA_HEADER_HOST
  void BuildLogPdfTable(int Z,
                        const double xmin,
                        const double xmax,
                        const int nrow,
                        const int ncol,
                        double *p);

public:
  // Auxiliary methods

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasSampler* GetSampler() {return fAliasSampler;}

  VECPHYS_CUDA_HEADER_BOTH
  void SetSampler(GUAliasSampler* sampler) { fAliasSampler = sampler ;}

private: 
  // Implementation methods 
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void 
  CrossSectionKernel(typename Backend::double  energyIn,
                     typename Backend::Index_t   zElement,
                     typename Backend::double& sgimaOut);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void 
  InteractKernel(typename Backend::double energyIn, 
                 typename Backend::Index_t   zElement,
                 typename Backend::double& energyOut,
                 typename Backend::double& sinTheta);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::double
  SampleSinTheta(typename Backend::double energyIn,
                 typename Backend::double energyOut) const; 

  VECPHYS_CUDA_HEADER_BOTH 
  void SampleByCompositionRejection(double energyIn,
                                    double& energyOut,
                                    double& sinTheta);

  VECPHYS_CUDA_HEADER_BOTH void 
  GetG4CrossSection(double  energyIn, 
                    const int zElement,
                    double& sgimaOut);
  
  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<ComptonKleinNishina>;

private:
  GUAliasSampler* fAliasSampler; 
  Precision fMinX;   // E Minimum - lowest energy for projectile
  Precision fMaxX;
  Precision fMaxZelement; // 
  
  //Sampling Tables
  int fNrow;
  int fNcol;
};

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
ComptonKleinNishina::CrossSectionKernel(typename Backend::double  energy, 
                                        typename Backend::Index_t   Z,
                                        typename Backend::double& sigmaOut)
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::double double;

  sigmaOut = 0.;
  Bool_t belowLimit = Bool_t(false);
  //low energy limit
  belowLimit |= ( energy < fLowEnergyLimit );
  if(Backend::early_returns && IsFull(belowLimit)) return;  

  double Z2 = Z*Z;
  double p1 =  2.7965e-1 +  1.9756e-5*Z + -3.9178e-7*Z2;
  double p2 = -1.8300e-1 + -1.0205e-2*Z +  6.8241e-5*Z2;
  double p3 =  6.7527    + -7.3913e-2*Z +  6.0480e-5*Z2;
  double p4 = -1.9798e+1 +  2.7079e-2*Z +  3.0274e-4*Z2;

  Bool_t condZ = (Z < 1.5);
  double T0 = 0.0; 
  CondAssign(condZ, 15.*keV, 40.*keV, &T0);  

  double X  =  Max(energy,T0)/electron_mass_c2;
  double X2 = X*X;
  double sigma = p1*Log(1.+2.*X)/X
          + (p2 + p3*X + p4*X2)/(1. + 20.*X + 230.*X2 + 440.*X2*X);
  sigmaOut = Z*sigma*barn;

  Bool_t condE = Bool_t(false);
  condE |= (energy > T0);
  if(Backend::early_returns && IsFull(condE)) return;  

  //correction when energy < T0
  double dT0 = 1.*keV;
  X = (T0+dT0) / electron_mass_c2 ;
  sigma = p1*log(1.+2.*X)/X
          + (p2 + p3*X + p4*X2)/(1. + 20.*X + 230.*X2 + 440.*X2*X);
  
  double   c1 = -T0*(Z*sigma*barn-sigmaOut)/(sigmaOut*dT0);
  double   c2 = 0.150;
  MaskedAssign( !condZ, 0.375-0.0556*Log(1.*Z) , &c2 );  
  double    y = Log(energy/T0);
  MaskedAssign(!condE, sigmaOut*Exp(-y*(c1+c2*y)),&sigmaOut);

  //this is the case if one of E < belowLimit 
  MaskedAssign(belowLimit, 0.0,&sigmaOut);
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
ComptonKleinNishina::InteractKernel(typename Backend::double  energyIn, 
                                    typename Backend::Index_t   zElement,
                                    typename Backend::double& energyOut,
                                    typename Backend::double& sinTheta)
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::double double;

  Index_t   index;
  Index_t   icol;
  double  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,index,icol,fraction);

  double probNA;
  double aliasInd;

  //this did not used to work - Fixed SW
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);
  
  double mininumE = energyIn/(1+2.0*energyIn*inv_electron_mass_c2);
  double deltaE = energyIn - mininumE;

  energyOut = mininumE + fAliasSampler->SampleX<Backend>(deltaE,probNA,
                                                aliasInd,icol,fraction);
  sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);
   
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::double 
ComptonKleinNishina::SampleSinTheta(typename Backend::double energyIn,
                                    typename Backend::double energyOut) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::double double;

  //angle of the scatterred photon

  double epsilon = energyOut/energyIn;

  Bool_t condition = epsilon > 1.0;

  MaskedAssign( condition, 1.0 , &epsilon );

  double E0_m    = inv_electron_mass_c2*energyIn;
  double onecost = (1.0 - epsilon)/(epsilon*E0_m);
  double sint2   = onecost*(2.-onecost);

  double sinTheta = 0.5;
  Bool_t condition2 = sint2 < 0.0;

  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );   
  
  return sinTheta;
}

} // end namespace impl
} // end namespace vecphys

#endif
