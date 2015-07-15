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
  void SampleByCompositionRejection(double energyIn,
                                    double& energyOut,
                                    double& sinTheta);
  
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
ComptonKleinNishina::InteractKernel(typename Backend::Double_t  energyIn, 
                                    typename Backend::Index_t   zElement,
                                    typename Backend::Double_t& energyOut,
                                    typename Backend::Double_t& sinTheta)
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
