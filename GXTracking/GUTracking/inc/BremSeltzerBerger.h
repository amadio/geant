#ifndef BremSeltzerBerger_H
#define BremSeltzerBerger_H 1

#include "backend/Backend.h"
#include "base/PhysicalConstants.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "Physics2DVector.h"
#include <fstream>

#include "EmModelBase.h"

namespace vecphys {

VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasSampler; )
inline namespace VECPHYS_IMPL_NAMESPACE {

class BremSeltzerBerger : public EmModelBase<BremSeltzerBerger>
{
public:

  VECPHYS_CUDA_HEADER_HOST
  BremSeltzerBerger(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  BremSeltzerBerger(Random_t* states, int threadId, 
                GUAliasSampler* sampler, Physics2DVector* sbData); 

  VECPHYS_CUDA_HEADER_BOTH
  ~BremSeltzerBerger();

  // Initializes this class and its sampler 
  VECPHYS_CUDA_HEADER_HOST
  void BuildOneTable(int Z,
                     const double xmin,
                     const double xmax,
                     const int nrow,
                     const int ncol);

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

  VECPHYS_CUDA_HEADER_BOTH
  Physics2DVector* GetSBData() {return fDataSB;}

  VECPHYS_CUDA_HEADER_HOST bool
  RetrieveSeltzerBergerData(std::ifstream& in, Physics2DVector *vec2D);

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
  SampleSinTheta(typename Backend::Double_t energyIn) const;

  VECPHYS_CUDA_HEADER_BOTH
  void SampleByCompositionRejection(int     elementZ,
                                    double  energyIn,
                                    double& energyOut,
                                    double& sinTheta);

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection(int Zelement, double Ein, double Eout) const;

  VECPHYS_CUDA_HEADER_BOTH double
  GetG4CrossSection(double  energyIn, 
                    const int zElement);

  VECPHYS_CUDA_HEADER_BOTH double
  ComputeXSectionPerAtom(double cut);

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<BremSeltzerBerger>;

private:
  GUAliasSampler* fAliasSampler; 
  Physics2DVector* fDataSB;
  Precision fMinX;
  Precision fMaxX;
  Precision fMaxZelement; 

  //Sampling Tables
  int fNrow;
  int fNcol;
};

//Implementation
template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
BremSeltzerBerger::CrossSectionKernel(typename Backend::Double_t  energy, 
                                      typename Backend::Index_t   Z)
{
  return 1.0;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
BremSeltzerBerger::InteractKernel(typename Backend::Double_t  energyIn,
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

  //Seltzer-Berger specific
  
  // To-do: apply densityFactor (dummy for now) and calculate deltaY correctly
  // densityFactor = (Migdal constant)x(electron density of the material); 
  Double_t densityFactor = 1.0;
  
  Double_t emin = Min(fMinX, energyIn);
  Double_t emax = Min(fMaxX, energyIn);

  Double_t totalEnergy = energyIn + electron_mass_c2;
  Double_t densityCorr = densityFactor*totalEnergy*totalEnergy;
  Double_t minY = Log(emin*emin + densityCorr);
  Double_t maxY = Log(emax*emax + densityCorr);
  Double_t deltaY = maxY - minY;

  Double_t yhat = fAliasSampler->SampleX<Backend>(deltaY,probNA,
                                                  aliasInd,icol,fraction);

  energyOut =  Sqrt(Max(Exp(minY + yhat)- densityCorr,0.0));
  sinTheta = SampleSinTheta<Backend>(energyOut);
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t 
BremSeltzerBerger::SampleSinTheta(typename Backend::Double_t energyIn) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  //angle of the radiated photon 
  //based on G4DipBustGenerator::SampleDirection
  
  Double_t c = 4. - 8.*UniformRandom<Backend>(fRandomState,fThreadId);
  Double_t a;
  Double_t signc; 
  Bool_t condition = c > 0.;
  MaskedAssign(  condition,  1. , &signc );
  MaskedAssign( !condition, -1. , &signc );
  MaskedAssign(  condition,  c , &a );
  MaskedAssign( !condition, -c , &a );

  Double_t delta  = Sqrt(a*a+4.);
  delta += a;
  delta *= 0.5; 

  //To-do:  Vc does not support pow 
  //  Double_t cofA = -signc*Pow(delta, 1./3.);
  Double_t cofA = -signc*Sqrt(delta); //temporary replace Sqrt by pow

  Double_t cosTheta = cofA - 1./cofA;

  Double_t tau  = energyIn/electron_mass_c2;
  Double_t beta = Sqrt(tau*(tau + 2.))/(tau + 1.);

  cosTheta = (cosTheta + beta)/(1 + cosTheta*beta);

  Double_t sinTheta = Sqrt((1 - cosTheta)*(1 + cosTheta));

  return sinTheta;
}

} // end namespace impl
} // end namespace vecphys

#endif