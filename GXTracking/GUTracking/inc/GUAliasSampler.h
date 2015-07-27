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
#include "GUAliasTableManager.h"

namespace vecphys {

  //VECPHYS_DEVICE_DECLARE_CONV( GUAliasSampler )

  //class GUAliasSampler;

inline namespace VECPHYS_IMPL_NAMESPACE
{

class GUAliasSampler
{
public: 

  VECPHYS_CUDA_HEADER_HOST
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
                 GUAliasTableManager* table 
                 );  
  
  VECPHYS_CUDA_HEADER_BOTH
  ~GUAliasSampler();

  VECPHYS_CUDA_HEADER_BOTH
  void PrintTable();

  int GetMaxZ() const { return fMaxZelement; }

  VECPHYS_CUDA_HEADER_HOST
  void BuildAliasTable( int z, int nrow, int ncol, const double *pdf );

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasTableManager* GetAliasTableManager(){ return fAliasTableManager ;}

  // Backend Implementation:
  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void SampleBin( typename Backend::double  kineticEnergy,
                  typename Backend::Index_t   &index,
                  typename Backend::Index_t   &icol,
                  typename Backend::double  &fraction) const;

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void SampleLogBin( typename Backend::double  kineticEnergy,
                     typename Backend::Index_t   &index,
                     typename Backend::Index_t   &icol,
                     typename Backend::double  &fraction) const;

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::double
  SampleX(typename Backend::double rangeSampled, 
          typename Backend::double probNA,   
          typename Backend::double aliasInd, 
          typename Backend::Index_t  icol,     
          typename Backend::double fraction  
          ) const;


  template<class Backend>
  inline
  VECPHYS_CUDA_HEADER_BOTH
  void
  GatherAlias(typename Backend::Index_t   index, 
              typename Backend::Index_t   zElement,
              typename Backend::double &probNA,  
              typename Backend::double &aliasInd 
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
  double   fInverseLogBinIncoming;
  
  //  For the sampled variable 
  const    int fSampledNumEntries;   //  Old name fNcol  (number of Columns)
  double   fInverseBinSampled; 
  GUAliasTableManager* fAliasTableManager; 
};

// Backend Implementation

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void GUAliasSampler::
SampleBin(typename Backend::double kineticEnergy,
          typename Backend::Index_t  &index,    // ~ sampled value
          typename Backend::Index_t  &icol,     // ~ input Energy
          typename Backend::double &fraction  //  in sampled variable
         ) const
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::double double;
  typedef typename Backend::Bool_t Bool_t;
  typedef typename Backend::int  int;

  //select the alias table for incoming energy 
  double eloc  = (kineticEnergy - fIncomingMin)*fInverseBinIncoming;
  Index_t  irow  = Floor(eloc);
  double efrac = eloc -1.0*irow;  
  double u1 = UniformRandom<Backend>(fRandomState,int(fThreadId));

  Bool_t condition = u1 <= efrac ;
  // if branch
  MaskedAssign( condition, irow , &irow ); // at the lower edge
  MaskedAssign( condition, irow + 1 , &irow ); // at the upper edge

  double r1 = (fSampledNumEntries-1)*UniformRandom<Backend>(fRandomState,int(fThreadId));
  icol = Floor(r1);
  fraction = r1 - 1.0*icol;

  // Was rangeSampled /(fSampledNumEntries-1);
  index = irow*fSampledNumEntries  + icol;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
void GUAliasSampler::
SampleLogBin(typename Backend::double kineticEnergy,
             typename Backend::Index_t  &index,    // ~ sampled value
             typename Backend::Index_t  &icol,     // ~ input Energy
             typename Backend::double &fraction  //  in sampled variable
             ) const
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::double double;
  typedef typename Backend::Bool_t Bool_t;
  typedef typename Backend::int  int;

  //select the alias table for incoming energy 
  double eloc = (Log(kineticEnergy) - Log(fIncomingMin))*fInverseLogBinIncoming;
  Index_t  irow = Floor(eloc);
  double efrac = eloc -1.0*irow;  
  
  double u1 = UniformRandom<Backend>(fRandomState,int(fThreadId));

  Bool_t condition = u1 <= efrac ;

  MaskedAssign( condition, irow , &irow );     // at the lower edge 
  MaskedAssign( condition, irow + 1 , &irow ); // at the upper edge

  //select the sampling bin
  double r1 = (fSampledNumEntries-1)*UniformRandom<Backend>(fRandomState,int(fThreadId));
  icol = Floor(r1);
  fraction = r1 - 1.0*icol;

  // Was rangeSampled /(fSampledNumEntries-1);
  index = irow*fSampledNumEntries  + icol;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::double
GUAliasSampler::
SampleX(typename Backend::double rangeSampled, 
        typename Backend::double probNA,   
        typename Backend::double aliasInd, 
        typename Backend::Index_t  icol,     
        typename Backend::double fraction  
       ) const
{
  typedef typename Backend::int    int;
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::double double;

  double r1 = UniformRandom<Backend>(fRandomState,int(fThreadId));

  Bool_t condition = r1 <= probNA;
  double xd, xu;
  double binSampled = rangeSampled * fInverseBinSampled;

  // if branch

  MaskedAssign( condition, icol*binSampled ,     &xd );   // Stores into xd
  MaskedAssign( condition, (icol+1)*binSampled , &xu );   //        into xu

  // else branch

  MaskedAssign( !condition,  aliasInd*binSampled    , &xd );
  MaskedAssign( !condition, (aliasInd+1)*binSampled , &xu );

  double x = (1 - fraction) * xd + fraction* xu;

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
            typename Backend::double  &probNA,
            typename Backend::double  &aliasInd
           ) const
{
#ifdef CHECK
  if( zElement <= 0  || zElement > fMaxZelement )
  {
    printf(" Illegal zElement = %d\n",zElement);
  }
#endif
  //  assert( (zElement > 0)  && (zElement <= fMaxZelement) );

  int     intIndex= (int) index;
   
#ifdef CHECK
  //  int     tableSize= fAliasTable[zElement]->SizeOfGrid();
  //  INDEX_CHECK( intIndex, 0, tableSize, "Index", "TableSize" );
#endif
  //  assert( (intIndex >= 0) && (intIndex < tableSize) );

  int tableIndex = fAliasTableManager->GetTableIndex(zElement);
  probNA=   (fAliasTableManager->GetAliasTable(tableIndex))->fProbQ[ intIndex ];
  aliasInd= (fAliasTableManager->GetAliasTable(tableIndex))->fAlias[ intIndex ];
}

// Specialisation for all vector-type backends - Vc for now
#ifndef VECPHYS_NVCC
template<>
inline
VECPHYS_CUDA_HEADER_BOTH
void GUAliasSampler::
GatherAlias<kVc>(typename kVc::Index_t    index, 
                 typename kVc::Index_t    zElement,
                 typename kVc::double  &probNA,  
                 typename kVc::double  &aliasInd
                ) const 
{
  //gather for alias table lookups - (backend type has no ptr arithmetic)
  for(int i = 0; i < kVc::kSize ; ++i) 
  {
    int z= zElement[i];
    int ind = index[i];

    if(ind < 0) {
      // printf("Warning: negative index! - in GUPhotoElectronSauterGavrila\n");
      ind = 0;
    }

    assert( z > 0  && z <= fMaxZelement );
    //    assert( ind >= 0 && ind < fAliasTable[z]->SizeOfGrid() );

    int tableIndex = fAliasTableManager->GetTableIndex(z);
    probNA[i]=   (fAliasTableManager->GetAliasTable(tableIndex))->fProbQ[ ind ];
    aliasInd[i]= (fAliasTableManager->GetAliasTable(tableIndex))->fAlias[ ind ];
  }
}
#endif

} // end namespace impl
} // end namespace vecphys

#endif
