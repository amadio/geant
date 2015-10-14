#ifndef ComptonKleinNishina_H
#define ComptonKleinNishina_H 1

#include "backend/Backend.h"
#include "base/PhysicalConstants.h"

#include "GUTrack.h"
#include "GUAliasSampler.h"

#include "EmModelBase.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class ComptonKleinNishina : public EmModelBase<ComptonKleinNishina>
{
public:

  VECPHYS_CUDA_HEADER_HOST
  ComptonKleinNishina(Random_t* states = 0, int threadId = -1);

  VECPHYS_CUDA_HEADER_BOTH
  ComptonKleinNishina(Random_t* states, int threadId, GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH 
  ~ComptonKleinNishina(){}

  //interfaces for tables
  VECPHYS_CUDA_HEADER_HOST 
  void BuildCrossSectionTablePerAtom(int Z);

  VECPHYS_CUDA_HEADER_HOST
  void BuildPdfTable(int Z, double *p);

public:
  // Auxiliary methods
  VECPHYS_CUDA_HEADER_BOTH
  GUAliasSampler* GetSampler() {return fAliasSampler;}

  VECPHYS_CUDA_HEADER_BOTH
  void SetSampler(GUAliasSampler* sampler) { fAliasSampler = sampler ;}


private: 
  // Implementation methods 
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH typename 
  Backend::Double_t
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

  VECPHYS_CUDA_HEADER_BOTH
  GUTrack& GetSecondaryElectron() { return fSecondaryElectron; }

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<ComptonKleinNishina>;

private:

  GUTrack fSecondaryElectron;
};

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
ComptonKleinNishina::CrossSectionKernel(typename Backend::Double_t  energy, 
                                        typename Backend::Index_t   Z)
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  Double_t sigmaOut = 0.;
  Bool_t belowLimit = Bool_t(false);
  //low energy limit
  belowLimit |= ( energy < fLowEnergyLimit );
  if(Backend::early_returns && IsFull(belowLimit)) return sigmaOut;  

  Double_t Z2 = Z*Z;
  Double_t p1 =  2.7965e-1 +  1.9756e-5*Z + -3.9178e-7*Z2;
  Double_t p2 = -1.8300e-1 + -1.0205e-2*Z +  6.8241e-5*Z2;
  Double_t p3 =  6.7527    + -7.3913e-2*Z +  6.0480e-5*Z2;
  Double_t p4 = -1.9798e+1 +  2.7079e-2*Z +  3.0274e-4*Z2;

  Bool_t condZ = (Z < 1.5);
  Double_t T0 = 0.0; 
  CondAssign(condZ, 15.*keV, 40.*keV, &T0);  

  Double_t X  =  Max(energy,T0)/electron_mass_c2;
  Double_t X2 = X*X;
  Double_t sigma = p1*Log(1.+2.*X)/X
          + (p2 + p3*X + p4*X2)/(1. + 20.*X + 230.*X2 + 440.*X2*X);
  sigmaOut = Z*sigma*barn;

  Bool_t condE = Bool_t(false);
  condE |= (energy > T0);
  if(Backend::early_returns && IsFull(condE)) return sigmaOut;  

  //correction when energy < T0
  Double_t dT0 = 1.*keV;
  X = (T0+dT0) / electron_mass_c2 ;
  sigma = p1*log(1.+2.*X)/X
          + (p2 + p3*X + p4*X2)/(1. + 20.*X + 230.*X2 + 440.*X2*X);
  
  Double_t   c1 = -T0*(Z*sigma*barn-sigmaOut)/(sigmaOut*dT0);
  Double_t   c2 = 0.150;
  MaskedAssign( !condZ, 0.375-0.0556*Log(1.*Z) , &c2 );  
  Double_t    y = Log(energy/T0);
  MaskedAssign(!condE, sigmaOut*Exp(-y*(c1+c2*y)),&sigmaOut);

  //this is the case if one of E < belowLimit 
  MaskedAssign(belowLimit, 0.0,&sigmaOut);
  return  sigmaOut;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
ComptonKleinNishina::InteractKernel(typename Backend::Double_t  energyIn, 
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
  
  Double_t mininumE = energyIn/(1+2.0*energyIn*inv_electron_mass_c2);
  Double_t deltaE = energyIn - mininumE;

  energyOut = mininumE + fAliasSampler->SampleXL<Backend>(zElement,
                                        deltaE,probNA,aliasInd,irow,icol);
  sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);

  //create the secondary electron

  //update the primary

  //  printf("icol = %d energyOut = %f %f %f %f\n",icol,energyOut,deltaE,aliasInd,probNA);   
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t 
ComptonKleinNishina::SampleSinTheta(typename Backend::Double_t energyIn,
                                    typename Backend::Double_t energyOut) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  //angle of the scatterred photon

  Double_t epsilon = energyOut/energyIn;

  Bool_t condition = epsilon > 1.0;

  MaskedAssign( condition, 1.0 , &epsilon );

  Double_t E0_m    = inv_electron_mass_c2*energyIn;
  Double_t onecost = (1.0 - epsilon)/(epsilon*E0_m);
  Double_t sint2   = onecost*(2.-onecost);

  Double_t sinTheta = 0.5;
  Bool_t condition2 = sint2 < 0.0;

  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );   
  
  return sinTheta;
}

} // end namespace impl
} // end namespace vecphys

#endif
