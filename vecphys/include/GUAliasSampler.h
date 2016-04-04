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

#include "GUAliasTable.h"
#include "GUAliasTableManager.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class GUAliasSampler
{
public:

  VECCORE_CUDA_HOST
  GUAliasSampler(Random_t* states,
                 int       threadId,
                 double    incomingMin,
                 double    incomingMax,
                 int       numEntriesIncoming,  // 'energy' (or log) of projectile
                 int       numEntriesSampled
                 );

  VECCORE_CUDA_HOST_DEVICE
  GUAliasSampler(Random_t* states, int threadId,
                 double incomingMin,
                 double incomingMax,
                 int    numEntriesIncoming, // 'energy' (or log) of projectile
                 int    numEntriesSampled,
                 GUAliasTableManager* table
                 );

  VECCORE_CUDA_HOST_DEVICE
  ~GUAliasSampler();

  VECCORE_CUDA_HOST_DEVICE
  void PrintTable();

  VECCORE_CUDA_HOST
  void BuildAliasTable( int z, const double *pdf );

  VECCORE_CUDA_HOST_DEVICE
  GUAliasTableManager* GetAliasTableManager(){ return fAliasTableManager ;}

  // Backend Implementation:
  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  void SampleBin( typename Backend::Double_v  kineticEnergy,
                  Index_v<typename Backend::Double_v>   &index,
                  Index_v<typename Backend::Double_v>   &icol,
                  typename Backend::Double_v  &fraction) const;

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  void SampleLogBin( typename Backend::Double_v  kineticEnergy,
                     Index_v<typename Backend::Double_v>   &irow,
                     Index_v<typename Backend::Double_v>   &icol,
                     typename Backend::Double_v  &fraction) const;

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  SampleX(typename Backend::Double_v rangeSampled,
          typename Backend::Double_v probNA,
          Index_v<typename Backend::Double_v> aliasInd,
          Index_v<typename Backend::Double_v>  icol,
          typename Backend::Double_v fraction ) const;

  template<class Backend>
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  SampleXL(Index_v<typename Backend::Double_v>  zElement,
           typename Backend::Double_v rangeSampled,
           typename Backend::Double_v probNA,
           Index_v<typename Backend::Double_v> aliasInd,
           Index_v<typename Backend::Double_v>  irow,
           Index_v<typename Backend::Double_v>  icol) const;

  template<class Backend>
  inline
  VECCORE_CUDA_HOST_DEVICE
  void
  GatherAlias(Index_v<typename Backend::Double_v>   index,
              Index_v<typename Backend::Double_v>   zElement,
              typename Backend::Double_v &probNA,
              Index_v<typename Backend::Double_v>  &aliasInd ) const;

  template<class Backend>
  inline
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  GetPDF(Index_v<typename Backend::Double_v> zElement,
         Index_v<typename Backend::Double_v> irow,
         Index_v<typename Backend::Double_v> icol ) const;

  //For atomic independent models
  template<class Backend>
  inline
  VECCORE_CUDA_HOST_DEVICE
  void
  GatherAlias(Index_v<typename Backend::Double_v>   index,
              typename Backend::Double_v &probNA,
              Index_v<typename Backend::Double_v>  &aliasInd ) const;

  template<class Backend>
  inline
  VECCORE_CUDA_HOST_DEVICE
  typename Backend::Double_v
  GetPDF(Index_v<typename Backend::Double_v> irow,
         Index_v<typename Backend::Double_v> icol ) const;

  //accessors
  VECCORE_CUDA_HOST_DEVICE
  double GetIncomingMin()  const { return fIncomingMin ; }

  VECCORE_CUDA_HOST_DEVICE
  double GetIncomingMax()  const { return fIncomingMax ; }

  VECCORE_CUDA_HOST_DEVICE
  int GetNumEntries()      const { return fInNumEntries; }

  VECCORE_CUDA_HOST_DEVICE
  int GetSamplesPerEntry() const { return fSampledNumEntries;}

private:
  Random_t* fRandomState;
  int       fThreadId;

  double   fIncomingMin; // Min of Incoming - e.g. e_Kinetic or Log(E_kinetic)
  double   fIncomingMax; // Max
  int      fInNumEntries;
  double   fLogIncomingMin;
  double   fInverseBinIncoming;
  double   fInverseLogBinIncoming;

  //  For the sampled variable
  const    int fSampledNumEntries;   //  Old name fNcol  (number of Columns)
  double   fInverseBinSampled;
  GUAliasTableManager* fAliasTableManager;
};

// Backend Implementation

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
void GUAliasSampler::
SampleBin(typename Backend::Double_v kineticEnergy,
          Index_v<typename Backend::Double_v>  &index,    // ~ sampled value
          Index_v<typename Backend::Double_v>  &icol,     // ~ input Energy
          typename Backend::Double_v &fraction  //  in sampled variable
         ) const
{
  typedef Index_v<typename Backend::Double_v>  Index_v<Double_v>;
  using Double_v = typename Backend::Double_v;
  typedef Mask_v<typename Backend::Double_v> Mask_v<Double_v>;
  typedef typename Backend::Int_t  Int_t;

  //select the alias table for incoming energy
  Double_v eloc  = (kineticEnergy - fIncomingMin)*fInverseBinIncoming;
  Index_v<Double_v>  irow  = Floor(eloc);
  Double_v efrac = eloc -1.0*irow;
  // to use fPower2Divisor
  //  Double_v eloc  = (kineticEnergy - fIncomingMin);
  //  Index_v<Double_v> irow = fPower2Divisor->GetBin<Backend>(eloc);
  //  Double_v efrac = fPower2Divisor->FractionWithinBin<Backend>(eloc,irow);

  Double_v u1 = UniformRandom<Backend>(fRandomState,Int_t(fThreadId));

  Mask_v<Double_v> useHigh = (u1 <= efrac) ;

  // irow = useHigh ? irow+1 : irow;
  MaskedAssign( useHigh, irow + 1 , &irow ); // at the upper edge

  Double_v r1 = fSampledNumEntries*UniformRandom<Backend>(fRandomState,Int_t(fThreadId));

  // Prepare output values
  icol = Floor(r1);
  fraction = r1 - 1.0*icol;

  // index = irow*fSampledNumEntries  + icol;  // No need to compute - no longer an output
}

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
void GUAliasSampler::
SampleLogBin(typename Backend::Double_v kineticEnergy,
             Index_v<typename Backend::Double_v>  &irow,     // input energy
             Index_v<typename Backend::Double_v>  &icol,     // sampled value
             typename Backend::Double_v &fraction  // within the sampled bin
             ) const
{
  using Double_v = typename Backend::Double_v;
  typedef Mask_v<typename Backend::Double_v> Mask_v<Double_v>;
  typedef typename Backend::Int_t  Int_t;

  //select the alias table for incoming energy
  Double_v eloc = (Log(kineticEnergy) - fLogIncomingMin)*fInverseLogBinIncoming;
  irow = Floor(eloc);
  Double_v efrac = eloc -1.0*irow;

  Double_v u1 = UniformRandom<Backend>(fRandomState,Int_t(fThreadId));

  Mask_v<Double_v> useHigh = (u1 <= efrac) ;

  // MaskedAssign( condition, irow , &irow );     // at the lower edge   --- Null operation!
  // MaskedAssign( condition, irow + 1 , &irow ); // at the upper edge

  // irow = useHigh ? irow+1 : irow;
  MaskedAssign( useHigh, irow + 1 , &irow ); // at the upper edge

  //select the sampling bin
  Double_v r1 = fSampledNumEntries*UniformRandom<Backend>(fRandomState,Int_t(fThreadId));
  icol = Floor(r1);
  fraction = r1 - 1.0*icol;

  // Was rangeSampled /(fSampledNumEntries-1);
  //  index = irow*fSampledNumEntries  + icol;
}

//    Sample distribution of secondary's 'X' - typically Energy
//      for given zElement ...
//    Feature of this method:  flat distribution within bin

template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
GUAliasSampler::
SampleX(typename Backend::Double_v rangeSampled,
        typename Backend::Double_v probNA,
        Index_v<typename Backend::Double_v>  aliasInd,
        Index_v<typename Backend::Double_v>  icol,
        typename Backend::Double_v fraction
       ) const
{
  typedef typename Backend::Int_t    Int_t;
  typedef Mask_v<typename Backend::Double_v>   Mask_v<Double_v>;
  using Double_v = typename Backend::Double_v;
  typedef Index_v<typename Backend::Double_v>  Index_v<Double_v>;

  Double_v r1 = UniformRandom<Backend>(fRandomState,Int_t(fThreadId));

  Mask_v<Double_v> useDirect = r1 <= probNA;  // Was Boot_t condition = ...
  Double_v xd, xu;
  Double_v binSampled = rangeSampled * fInverseBinSampled;

  Index_v<Double_v>       icolDist= icol;
  MaskedAssign( !useDirect, aliasInd, &icolDist );

  xd = icolDist*binSampled;
  xu = xd + binSampled;

  Double_v x = (1 - fraction) * xd + fraction* xu;

  return x;
}

//
//  Lacks a description of what the method does.
//  Since it is a complex method, it will benefit significantly from it.

//  Draft description:
//    Sample distribution of secondary's 'X' - typically Energy
//      for given zElement ...
//    Feature of this method:  linear interpolation using 'PDF'
template<class Backend>
VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v
GUAliasSampler::
SampleXL(Index_v<typename Backend::Double_v>  zElement,
         typename Backend::Double_v rangeSampled,
         typename Backend::Double_v probNA,
         Index_v<typename Backend::Double_v>  aliasInd,
         Index_v<typename Backend::Double_v>  irow,
         Index_v<typename Backend::Double_v>  icol) const
{
  typedef typename Backend::Int_t    Int_t;
  typedef Mask_v<typename Backend::Double_v>   Mask_v<Double_v>;
  using Double_v = typename Backend::Double_v;
  typedef Index_v<typename Backend::Double_v>  Index_v<Double_v>;

  Double_v r1 = UniformRandom<Backend>(fRandomState,Int_t(fThreadId));

  Mask_v<Double_v> condition = r1 <= probNA;
  Double_v xd, xu;
  Double_v binSampled = rangeSampled * fInverseBinSampled;

  Index_v<Double_v>       icolDist= icol;
  MaskedAssign( !condition, aliasInd, &icolDist );
  xd = icolDist * binSampled;
  xu = xd + binSampled;

  // Flat distribution was
  //  Double_v x = (1 - fraction) * xd + fraction* xu;

  //Using pdf of linear interpolation within the sampling bin based on the pdf
  //linear interpolation within the sampling bin based on the pdf
  Double_v x(0.);
  Double_v pd(0.);
  Double_v pu(0.);

  pd = GetPDF<Backend>(irow,icolDist);
  pu = GetPDF<Backend>(irow,icolDist+1);

  //* Obtain 'x' in interval [xd, xu] using pdf from linear interpolation
  //    (x,y) from (xd, pd) and (xu, pu)
  //  - Uses two random numbers in order to avoid square root
  Double_v r2 = UniformRandom<Backend>(fRandomState,Int_t(fThreadId));
  Double_v r3 = UniformRandom<Backend>(fRandomState,Int_t(fThreadId));

  Mask_v<Double_v> below = r2*(pd+pu) < (1.-r3)*pd + r3*pu;;

  MaskedAssign(  below, (1.-r3)*xd + r3*xu , &x);
  MaskedAssign( !below, r3*xd + (1.-r3)*xu , &x);

  //- Can simplify to
  // MaskedAssign( !below, 1.0-r3, r3);      //  if (!below) { r3 = 1.0 - r3; }
  // x = (1.-r3)*xd + r3*xu;

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
GatherAlias(Index_v<typename Backend::Double_v>    index,
            Index_v<typename Backend::Double_v>    zElement,
            typename Backend::Double_v  &probNA,
            Index_v<typename Backend::Double_v>   &aliasInd
           ) const
{
#ifdef CHECK
  //if( zElement <= 0  || zElement > fMaxZelement )
  //{
  //  printf(" Illegal zElement = %d\n",zElement);
  //}
#endif
  //  assert( (zElement > 0)  && (zElement <= fMaxZelement) );

  int     intIndex= (int) index;

#ifdef CHECK
  //  int     tableSize= fAliasTable[zElement]->SizeOfGrid();
  //  INDEX_CHECK( intIndex, 0, tableSize, "Index", "TableSize" );
#endif
  //  assert( (intIndex >= 0) && (intIndex < tableSize) );

  probNA=   (fAliasTableManager->GetAliasTable(zElement))->fProbQ[ intIndex ];
  aliasInd= (fAliasTableManager->GetAliasTable(zElement))->fAlias[ intIndex ];
}

template<class Backend>
inline
typename Backend::Double_v
GUAliasSampler::GetPDF(Index_v<typename Backend::Double_v> zElement,
                       Index_v<typename Backend::Double_v> irow,
                       Index_v<typename Backend::Double_v> icol) const
{
  using Double_v = typename Backend::Double_v;

  int     intIndex= (int) (fSampledNumEntries*irow + icol);
  Double_v pdf= (fAliasTableManager->GetAliasTable(zElement))->fpdf[intIndex];

  return pdf;
}

//For atomic independent models

template<class Backend>
inline
void GUAliasSampler::
GatherAlias(Index_v<typename Backend::Double_v>    index,
            typename Backend::Double_v  &probNA,
            Index_v<typename Backend::Double_v>   &aliasInd
           ) const
{
  int     intIndex= (int) index;
  probNA =   (fAliasTableManager->GetAliasTable(0))->fProbQ[ intIndex ];
  aliasInd = (fAliasTableManager->GetAliasTable(0))->fAlias[ intIndex ];
}


template<class Backend>
inline
typename Backend::Double_v
GUAliasSampler::GetPDF(Index_v<typename Backend::Double_v> irow,
                       Index_v<typename Backend::Double_v> icol) const
{
  using Double_v = typename Backend::Double_v;

  int     intIndex= (int) (fSampledNumEntries*irow + icol);
  Double_v pdf= (fAliasTableManager->GetAliasTable(0))->fpdf[intIndex];

  return pdf;
}

// Specialisation for all vector-type backends - Vc for now
#ifndef VECCORE_NVCC
template<>
inline
VECCORE_CUDA_HOST_DEVICE
void GUAliasSampler::
GatherAlias<kVc>(typename kVc::Index_v<Double_v>    index,
                 typename kVc::Index_v<Double_v>    zElement,
                 typename kVc::Double_v  &probNA,
                 typename kVc::Index_v<Double_v>   &aliasInd
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

    //assert( z > 0  && z <= fMaxZelement );
    //    assert( ind >= 0 && ind < fAliasTable[z]->SizeOfGrid() );

    probNA[i]=   (fAliasTableManager->GetAliasTable(z))->fProbQ[ ind ];
    aliasInd[i]= (fAliasTableManager->GetAliasTable(z))->fAlias[ ind ];
  }
}

template<>
inline
VECCORE_CUDA_HOST_DEVICE
typename kVc::Double_v
GUAliasSampler::GetPDF<kVc>(typename kVc::Index_v<Double_v> zElement,
                            typename kVc::Index_v<Double_v> irow,
                            typename kVc::Index_v<Double_v> icol) const
{
  typedef typename kVc::Double_v Double_v;

  Double_v pdf;

  for(int i = 0; i < kVc::kSize ; ++i)
  {
    int z= zElement[i];
    int ind = fSampledNumEntries*irow[i] + icol[i];

    if(ind < 0) {
      ind = 0;
    }

    //assert( z > 0  && z <= fMaxZelement );
    pdf[i] = (fAliasTableManager->GetAliasTable(z))->fpdf[ ind ];
  }
  return pdf;
}

//For atomic indedepend models

template<>
inline
VECCORE_CUDA_HOST_DEVICE
void GUAliasSampler::
GatherAlias<kVc>(typename kVc::Index_v<Double_v>    index,
                 typename kVc::Double_v  &probNA,
                 typename kVc::Index_v<Double_v>   &aliasInd
                ) const
{
  //gather for alias table lookups - (backend type has no ptr arithmetic)
  for(int i = 0; i < kVc::kSize ; ++i)
  {
    int ind = index[i];
    probNA[i]=   (fAliasTableManager->GetAliasTable(0))->fProbQ[ ind ];
    aliasInd[i]= (fAliasTableManager->GetAliasTable(0))->fAlias[ ind ];
  }
}

template<>
inline
VECCORE_CUDA_HOST_DEVICE
typename kVc::Double_v
GUAliasSampler::GetPDF<kVc>(typename kVc::Index_v<Double_v> irow,
                            typename kVc::Index_v<Double_v> icol) const
{
  typedef typename kVc::Double_v Double_v;

  Double_v pdf;

  for(int i = 0; i < kVc::kSize ; ++i)
  {
    int ind = fSampledNumEntries*irow[i] + icol[i];
    pdf[i] = (fAliasTableManager->GetAliasTable(0))->fpdf[ ind ];
  }
  return pdf;
}


#endif

} // end namespace impl
} // end namespace vecphys

#endif
