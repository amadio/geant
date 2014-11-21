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
#include "backend/vc/Backend.h"

#include "GURandom.h"

using namespace VECGEOM_NAMESPACE;

class GUAliasSampler
{
public: 
  GUAliasSampler(int    Zelement, 
                 double incomingMin, 
                 double incomingMax,
                 int    numEntriesIncoming,  // for 'energy' (or log) of projectile
                 int    numEntriesSampled 
                 );  

  ~GUAliasSampler();
  void PrintTable();

  template<class Backend>
  typename Backend::Double_t
  Sample( typename Backend::Double_t energyIn, 
          typename Backend::Double_t deltaY    ) const;  

  /* Builds (part) of our tables ( for one ene)
   *
  */
   // Builds all our table - must be called during initialisation 
  //  void SetPdf( double const * pdf );
  //  void BuildRow( int numVal, double* differentialXSection, double xMin, double xMax );

  //private:
// Implementation methods: 
  void BuildAliasTables( const int nrow, const int ncol, double   *pdf );

  template<class Backend>
  typename Backend::Double_t
  SampleX( typename Backend::Index_t  irow,   // ~ sampled value
           typename Backend::Index_t  icol,   // ~ input Energy
           typename Backend::Double_t binSize,
           typename Backend::Double_t remainderX  //  in sampled variable
          ) const;

  template<class Backend>
  FQUALIFIER void
  GetBin( typename Backend::Double_t  kineticEnergy,
	  typename Backend::Index_t     &irow,
	  typename Backend::Index_t     &icol,
          typename Backend::Double_t  &t) const;

private:
  const int      fZelement; 
  
  const double   fIncomingMin; // Min of Incoming - e.g. e_Kinetic or Log(E_kinetic)
  const double   fIncomingMax; // Max
  const int      fInNumEntries;
  const double   fInverseBinIncoming;
  // const double   fInBinSize; 
  
  //  For the sampled variable 
  const int      fSampledNumEntries;   //  Old name fNcol  (number of Columns)
  const double   fInverseBinSampled; 
  const double   fSampledBinSize; 
  // Values needed per bin 
  //  -- Initialised how ??
  double*  fSampledMin; // Minimum value of 'x' sampled variable
  double*  fSampledMax; // Maximum value of 'x' sampled variable
  
  // Effective 2-dimensional arrays - size is fInNumEntries * fSampledNumEntries
  double * fpdf; // Original distribution
  
  double * fProbQ; // Non-alias probability ( sometimes called q in other code )
  int    * fAlias; // Alias table           -- could become Index_t ?
};

// implementation of template functions

template<class Backend>
typename Backend::Double_t
GUAliasSampler::
Sample( typename Backend::Double_t energyIn,
    typename Backend::Double_t deltaY ) const
{
  typedef typename Backend::Double_t Double_t;
  typedef typename Backend::Index_t  Index_t;

  Index_t   irow, icol;
  Double_t  fraction;

  GetBin<Backend>(energyIn, irow, icol, fraction);
  Double_t x = SampleX<Backend>( irow, icol, deltaY, fraction);

  return x;
}

template<class Backend>
FQUALIFIER void
GUAliasSampler::
GetBin(typename Backend::Double_t  kineticEnergy,
       typename Backend::Index_t   &irow,
       typename Backend::Index_t   &icol,
       typename Backend::Double_t  &t) const
{
  typedef typename Backend::Double_t Double_t;

  irow = VECGEOM_NAMESPACE::Floor((kineticEnergy - fIncomingMin)*fInverseBinIncoming);
  Double_t r1 = (fSampledNumEntries-1)*Double_t::Random();
  icol = VECGEOM_NAMESPACE::Floor(r1);
  t = r1 - 1.0*icol;
}

template<class Backend>
typename Backend::Double_t
GUAliasSampler::
SampleX( typename Backend::Index_t  irow,   // ~ sampled value
         typename Backend::Index_t  icol,   // ~ input Energy
         typename Backend::Double_t rangeSampled,
         typename Backend::Double_t remainderX  //  in sampled variable
        ) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

  Double_t r1 = Double_t::Random();//GUUniformRand(0,-1);

  Double_t xd, xu;
  Double_t binSampled = rangeSampled * fInverseBinSampled;
        // Was rangeSampled /(fSampledNumEntries-1);

  Double_t probNA;   // Non-alias probability
  Double_t aliasInd; //  This is really an integer -- could be Index_t !?
  // fill probNA from table
  //  Index_t index = irow*fSampledNumEntries  + icol;
  Index_t index = irow*fSampledNumEntries  + icol;

  // Gather
  for( int i=0;i < Double_t::Size; ++i )
  {
    //    int iEntry= ConverfractionoInteger(index[i]);
    //    probNA[i]=    fPDFY[ iEntry ]; // index[i] ];
    //    aliasInd[i]=  fPDFA[ iEntry ]; // index[i] ];
    
    probNA[i]=    fProbQ[ (int) index[i] ];
    aliasInd[i]=  fAlias[ (int) index[i] ];
  }
  // should investigate here whether gather is supported in Vc

  Bool_t condition = r1 <= probNA;

  // if branch

  MaskedAssign( condition, icol*binSampled ,     &xd );   // Stores into xd
  MaskedAssign( condition, (icol+1)*binSampled , &xu );   //        into xu

  // else branch

  MaskedAssign( !condition,  aliasInd*binSampled    , &xd );
  MaskedAssign( !condition, (aliasInd+1)*binSampled , &xu );

  Double_t x = (1 - remainderX) * xd + remainderX* xu;

  return x;
}

#endif
