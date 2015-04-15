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

#include "GUAliasTable.h"

#define   CHECK   1
//  Uncomment to provide checking of indices and information about their failures

namespace vecphys {

VECPHYS_DEVICE_DECLARE_CONV( GUAliasSampler )

  //VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasTable; )

inline namespace VECPHYS_IMPL_NAMESPACE
{

class GUAliasSampler
{
public: 

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasSampler(Random_t* states,
                 int       threadId,
                 int       maxZelement,   //  ==> Now For all Z
                 double    incomingMin, 
                 double    incomingMax,
                 int       numEntriesIncoming,  // 'energy' (or log) of projectile
                 int       numEntriesSampled 
                 );

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasSampler(Random_t* states, int threadId,
                 // int    Zelement,   //  ==> Now For all Z
                 double incomingMin, 
                 double incomingMax,
                 int    numEntriesIncoming, // 'energy' (or log) of projectile
                 int    numEntriesSampled, 
                 GUAliasTable* table 
                 );  
  
  VECPHYS_CUDA_HEADER_BOTH
  ~GUAliasSampler();

  VECPHYS_CUDA_HEADER_BOTH
  void PrintTable();

  VECPHYS_CUDA_HEADER_BOTH
  void GetAlias(int     index,
                int     zElement,
                double &probNA,  
                int    &aliasInd) const;

  int GetMaxZ() const { return fMaxZelement; }

  VECPHYS_CUDA_HEADER_BOTH
  void BuildAliasTable( int z, int nrow, int ncol, const double *pdf );

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasTable* GetAliasTable(int z){ return (fAliasTable && z < fMaxZelement && z>0) ? fAliasTable[z]: NULL ;}

  VECPHYS_CUDA_HEADER_BOTH
  void SetAliasTable(int z, GUAliasTable* table) { if( table && z < fMaxZelement && z>0) { delete fAliasTable[z]; fAliasTable[z] = table ;} }
  
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
  inline
  VECPHYS_CUDA_HEADER_BOTH
  void
  GatherAlias(typename Backend::Index_t   index, 
              typename Backend::Index_t   zElement,
              typename Backend::Double_t &probNA,  
              typename Backend::Double_t &aliasInd 
              ) const;

  int GetNumEntries()      const { return fInNumEntries; }      // 'Input' values:  E', log(E')
  int GetSamplesPerEntry() const { return fSampledNumEntries;}
  
private:
  Random_t*      fRandomState;
  int            fThreadId;
 
  int      fMaxZelement; 
  
  double   fIncomingMin; // Min of Incoming - e.g. e_Kinetic or Log(E_kinetic)
  double   fIncomingMax; // Max
  int      fInNumEntries;
  double   fInverseBinIncoming;
  
  //  For the sampled variable 
  const    int fSampledNumEntries;   //  Old name fNcol  (number of Columns)
  double   fInverseBinSampled; 
  // double   fSampledBinSize; 

  // double*  fSampledMin; // Minimum value of 'x' sampled variable
  // double*  fSampledMax; // Maximum value of 'x' sampled variable
  
  // Effective 2-dimensional arrays - size is fInNumEntries * fSampledNumEntries
  GUAliasTable** fAliasTable;    // GUAliasTable* fAliasTable[];
  // std::vector <GUAliasTable*>fAliasTable; 
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
  Double_t r1 = (fSampledNumEntries-1)*UniformRandom<Backend>(fRandomState,fThreadId);
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

  Double_t r1 = UniformRandom<Backend>(fRandomState,fThreadId);

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

#define INDEX_CHECK( AnIndex, BminVal, CmaxVal, DindexName, EmaxName ) \
   if( (AnIndex < BminVal) || (AnIndex > CmaxVal) ){                   \
     printf(" Illegal %s = %d vs min = %d and max ( %s ) = %d\n",      \
	    DindexName,AnIndex,BminVal,EmaxName,CmaxVal);              \
   }

// Scalar method - to be used ONLY for 'scalar-type' backends
//                  i.e. currently: scalar & CUDA
template<class Backend>
inline
void GUAliasSampler::
GatherAlias(typename Backend::Index_t    index,
            typename Backend::Index_t    zElement,
            typename Backend::Double_t  &probNA,
            typename Backend::Double_t  &aliasInd
           ) const
{
#ifdef CHECK
  if( zElement <= 0  || zElement > fMaxZelement )
  {
    printf(" Illegal zElement = %d\n",zElement);
  }
#endif
  assert( (zElement > 0)  && (zElement <= fMaxZelement) );
   
  int     intIndex= (int) index;
  int     tableSize= fAliasTable[zElement]->SizeOfGrid();
#ifdef CHECK
  // if( (intIndex < 0) || (intIndex > tableSize) )
  // {
  //    std::cout << " Illegal index = " << intIndex << " vs tablesize = " << tableSize << std::endl;
  // }
  
  //           index    <min  >max       name    name Max
  INDEX_CHECK( intIndex, 0, tableSize, "Index", "TableSize" );
#endif
  assert( (intIndex >= 0) && (intIndex < tableSize) );

  probNA=    fAliasTable[(int)zElement]->fProbQ[ intIndex ];
  aliasInd=  fAliasTable[(int)zElement]->fAlias[ intIndex ];
}

// Specialisation for all vector-type backends - Vc for now
#ifndef VECPHYS_NVCC
template<>
inline
VECPHYS_CUDA_HEADER_BOTH
void GUAliasSampler::
GatherAlias<kVc>(typename kVc::Index_t    index, 
                 typename kVc::Index_t    zElement,
                 typename kVc::Double_t  &probNA,  
                 typename kVc::Double_t  &aliasInd
                ) const 
{
  //gather for alias table lookups - (backend type has no ptr arithmetic)
  for(int i = 0; i < kVc::kSize ; ++i) 
  {
    int z= zElement[i];
    int ind = index[i];
    assert( z > 0  && z <= fMaxZelement );
    assert( ind >= 0 && ind < fAliasTable[z]->SizeOfGrid() );

    probNA[i]=    fAliasTable[z]->fProbQ[ ind ]; // (int) index[i] ];
    aliasInd[i]=  fAliasTable[z]->fAlias[ ind ]; // (int) index[i] ];
  }
}
#endif

} // end namespace impl
} // end namespace vecphys

#endif
