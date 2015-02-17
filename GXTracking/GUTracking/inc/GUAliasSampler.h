#ifndef GUAliasSampler_H
#define GUAliasSampler_H 1
// 
//  Alias Sampler template function implementation
//   For use in implementing methods for scalar, vector, GPU
//   Depends on Backend 'technology' of VecGeom
//
//  First version using 'Backend' - 22 Oct 2014
//   Authors:  Sandro Wenzel, Soon Y. Jun, John Apostolakis
//
//  Desing choice: One Sampler per element 
//                    (potentially later per composit material?)
//
//  First version assumes linear 'x' integrand
//   TODO: a measure will be required for the integration for log, theta etc.
// 
//  Based on first alias sampler by Soon Y. Jun - July 2014
#include "backend/Backend.h"

namespace vecphys {

inline namespace VECPHYS_IMPL_NAMESPACE {

class GUAliasSampler
{
public: 

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasSampler(Random_t* states, int threadId,
                 int    Zelement, 
                 double incomingMin, 
                 double incomingMax,
                 int    numEntriesIncoming, // 'energy' (or log) of projectile
                 int    numEntriesSampled 
                 );  

  VECPHYS_CUDA_HEADER_BOTH
  ~GUAliasSampler();

  VECPHYS_CUDA_HEADER_BOTH
  void PrintTable();

  VECPHYS_CUDA_HEADER_BOTH
  void GetAlias(int     index, 
                double &probNA,  
                int    &aliasInd) const;

  VECPHYS_CUDA_HEADER_BOTH
  void BuildAliasTables( const int nrow, const int ncol, double   *pdf );

  // Backend Implementation:
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void SampleBin( typename Backend::Double_t  kineticEnergy,
                  typename Backend::Index_t   &index,
                  typename Backend::Index_t   &icol,
                  typename Backend::Double_t  &fraction) const;

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  SampleX(typename Backend::Double_t rangeSampled, 
          typename Backend::Double_t probNA,   
          typename Backend::Double_t aliasInd, 
          typename Backend::Index_t  icol,     
          typename Backend::Double_t fraction  
          ) const;


  //#ifndef VECPHYS_NVCC
  template<class Backend>
  void
  GatherAlias(typename Backend::Index_t  index, 
              typename Backend::Double_t &probNA,  
              typename Backend::Double_t &aliasInd 
              ) const;

private:
  Random_t*      fRandomState;
  int            fThreadId;

  const int      fZelement; 
  
  const double   fIncomingMin; // Min of Incoming - e.g. e_Kinetic or Log(E_kinetic)
  const double   fIncomingMax; // Max
  const int      fInNumEntries;
  const double   fInverseBinIncoming;
  
  //  For the sampled variable 
  const int      fSampledNumEntries;   //  Old name fNcol  (number of Columns)
  const double   fInverseBinSampled; 
  const double   fSampledBinSize; 

  double*  fSampledMin; // Minimum value of 'x' sampled variable
  double*  fSampledMax; // Maximum value of 'x' sampled variable
  
  // Effective 2-dimensional arrays - size is fInNumEntries * fSampledNumEntries
  double * fpdf; // Original distribution
  
  double * fProbQ; // Non-alias probability ( sometimes called q in other code )
  int    * fAlias; // Alias table           -- could become Index_t ?
};

// Backend Implementation

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void GUAliasSampler::
SampleBin(typename Backend::Double_t kineticEnergy,
          typename Backend::Index_t  &index,    // ~ sampled value
          typename Backend::Index_t  &icol,     // ~ input Energy
          typename Backend::Double_t &fraction  //  in sampled variable
         ) const
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

  Index_t irow = Floor((kineticEnergy - fIncomingMin)*fInverseBinIncoming);
  Double_t r1 = (fSampledNumEntries-1)*UniformRandom(fRandomState,fThreadId);
  icol = Floor(r1);
  fraction = r1 - 1.0*icol;

  // Was rangeSampled /(fSampledNumEntries-1);
  index = irow*fSampledNumEntries  + icol;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
GUAliasSampler::
SampleX(typename Backend::Double_t rangeSampled, 
        typename Backend::Double_t probNA,   
        typename Backend::Double_t aliasInd, 
        typename Backend::Index_t  icol,     
        typename Backend::Double_t fraction  
       ) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;


  Double_t r1 = UniformRandom(fRandomState,fThreadId);

  Bool_t condition = r1 <= probNA;
  Double_t xd, xu;
  Double_t binSampled = rangeSampled * fInverseBinSampled;

  // if branch

  MaskedAssign( condition, icol*binSampled ,     &xd );   // Stores into xd
  MaskedAssign( condition, (icol+1)*binSampled , &xu );   //        into xu

  // else branch

  MaskedAssign( !condition,  aliasInd*binSampled    , &xd );
  MaskedAssign( !condition, (aliasInd+1)*binSampled , &xu );

  Double_t x = (1 - fraction) * xd + fraction* xu;

  return x;
}

#ifndef VECPHYS_NVCC

template<class Backend>
void GUAliasSampler::
GatherAlias(typename Backend::Index_t  index, 
            typename Backend::Double_t &probNA,  
            typename Backend::Double_t &aliasInd 
           ) const 
{
  //gather for alias table lookups - (backend type has no ptr arithmetic)
  for(int i = 0; i < Backend::kSize ; ++i) {
    probNA[i]=    fProbQ[ (int) index[i] ];
    aliasInd[i]=  fAlias[ (int) index[i] ];
  }
}
#endif

} // end namespace impl
} // end namespace vecphys

#endif
