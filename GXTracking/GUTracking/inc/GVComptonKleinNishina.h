#ifndef GVComptonKleinNishina_H
#define GVComptonKleinNishina_H 1

#include "backend/Backend.h"

#include "PhysicalConstants.h"
#include "GUTrack.h"

#include "GUAliasSampler.h"

#include "GUEmModelBase.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class GVComptonKleinNishina : public GUEmModelBase
{
public:

  VECPHYS_CUDA_HEADER_HOST
  GVComptonKleinNishina(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  GVComptonKleinNishina(Random_t* states, int threadId, 
                        GUAliasSampler* sampler); 

  VECPHYS_CUDA_HEADER_BOTH
  ~GVComptonKleinNishina();

  // Core method(s)
  // -------------------------------------------  
  VECPHYS_CUDA_HEADER_BOTH
  virtual void VInteractKernel(double energyIn, 
                               int   zElement,
                               double& energyOut,
                               double& sinTheta);

#ifndef VECPHYS_NVCC
  VECPHYS_CUDA_HEADER_BOTH
  virtual void VInteractKernel(kVc::Double_t  energyIn, 
			       kVc::Index_t  zElement,
			       kVc::Double_t& energyOut,
		               kVc::Double_t& sinTheta);
#endif

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void 
  InteractKernel(typename Backend::Double_t energyIn, 
                 typename Backend::Index_t zElement,
                 typename Backend::Double_t& energyOut,
                 typename Backend::Double_t& sinTheta) const;


  // Alternative Implementation method(s) - for reference/comparisons
  // ----------------------------------------------------------------
  VECPHYS_CUDA_HEADER_BOTH
  virtual void 
  SampleByCompositionRejection(double energyIn,
                               double& energyOut,
                               double& sinTheta);

  //  Initialisation methods
  // -------------------------------------------

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

public:

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  SampleSinTheta(typename Backend::Double_t energyIn,
                 typename Backend::Double_t energyOut) const; 

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t 
  TotalCrossSection(typename Backend::Double_t energy,
                    typename Backend::Double_t Z) const;

private: 
  // Implementation methods 

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

private:
  GUAliasSampler* fAliasSampler; 

  // Helper data members for GPU random -- to be replaced by use of a GPU manager class
  Random_t* fRandomState;
  int       fThreadId;

  Precision fMinX;   // E Minimum - lowest energy for projectile
  Precision fMaxX;
  // Precision fDeltaX;

  // Precision fMinY, fMaxY, fDeltaY;  // Energy limits for outgoing particles ? Not used

  Precision fMaxZelement; // 
  
  //Sampling Tables
  int fNrow;
  int fNcol;
};

//Implementation

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
GVComptonKleinNishina::InteractKernel(typename Backend::Double_t  energyIn, 
                                      typename Backend::Index_t   zElement,
                                      typename Backend::Double_t& energyOut,
                                      typename Backend::Double_t& sinTheta)
                                      const
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

  Index_t   index;
  Index_t   icol;
  Double_t  fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn,index,icol,fraction);

  Double_t probNA;
  Double_t aliasInd;

  //this did not used to work - Fixed SW
  fAliasSampler->GatherAlias<Backend>(index,zElement,probNA,aliasInd);
  
  Double_t mininumE = energyIn/(1+2.0*energyIn*inv_electron_mass_c2);
  Double_t deltaE = energyIn - mininumE;

  energyOut = mininumE + fAliasSampler->SampleX<Backend>(deltaE,probNA,
					        aliasInd,icol,fraction);
  sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t 
GVComptonKleinNishina::SampleSinTheta(typename Backend::Double_t energyIn,
                                      typename Backend::Double_t energyOut) 
                                      const
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

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t 
GVComptonKleinNishina::
TotalCrossSection(typename Backend::Double_t energy,
                  typename Backend::Double_t Z) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  Double_t sigma = 0.;

  //low energy limit
  Bool_t condition = energy > 10*keV;
  MaskedAssign( !condition, 0.0 , &sigma );

  if(Any(condition)) {
    Double_t Z2 =  Z*Z;
    //put coff's to a constant header
    Double_t p1 =  2.7965e-1 +  1.9756e-5*Z + -3.9178e-7*Z2;
    Double_t p2 = -1.8300e-1 + -1.0205e-2*Z +  6.8241e-5*Z2;
    Double_t p3 =  6.7527    + -7.3913e-2*Z +  6.0480e-5*Z2;
    Double_t p4 = -1.9798e+1 +  2.7079e-2*Z +  3.0274e-4*Z2;
 
    Double_t X = energy/electron_mass_c2;
    Double_t X2 = X*Z;

    Double_t tmpsigma = p1*log(1.+2.*X)/X
          + (p2 + p3*X + p4*X2)/(1. + 20.*X + 230.*X2 + 440.*X2*X);
    tmpsigma *= Z*barn;

    MaskedAssign( condition, tmpsigma , &sigma );
  }

  // 5% level improvements for low energy at 10-40keV for Hydrogen
  // and at 10-15keV for all Z    

  return sigma;
}

} // end namespace impl
} // end namespace vecphys

#endif
