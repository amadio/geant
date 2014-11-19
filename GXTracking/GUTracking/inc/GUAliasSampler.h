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
  template<class Backend>
  GUAliasSampler(int    Zelement, 
                 double incomingMin, 
                 double incomingMax,
                 int    numEntriesIncoming,  // for 'energy' (or log) of projectile
                 int    numEntriesSampled 
                 );  

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
  SampleX( typename Backend::Int_t  irow,   // ~ sampled value 
           typename Backend::Int_t  icol,   // ~ input Energy
           typename Backend::Double_t binSize,
           typename Backend::Double_t remainderX  //  in sampled variable
          );

  template<class Backend>
  FQUALIFIER void
  GetBin( typename Backend::Double_t  kineticEnergy,
	  typename Backend::Int_t     &irow, 
	  typename Backend::Int_t     &icol,
          typename Backend::Double_t  &t);

private:
  const int      fZelement; 
  
  const double   fIncomingMin; // Min of Incoming - e.g. e_Kinetic or Log(E_kinetic)
  const double   fIncomingMax; // Max
  const int      fInNumEntries;
  const double   fInverseBinIncoming;
  // const double   fInBinSize; 
  
  //  For the sampled variable 
  const int      fSampledNumEntries;   //  Old name fNcol  (number of Columns)
  const double   fSampledBinSize; 
  const double   fInverseBinSampled; 
  // Values needed per bin 
  //  -- Initialised how ??
  double*  fSampledMin; // Minimum value of 'x' sampled variable
  double*  fSampledMax; // Maximum value of 'x' sampled variable
  
  // Effective 2-dimensional arrays - size is fInNumEntries * fSampledNumEntries
  double * fpdf; // Original distribution
  
  double * fProbQ; // Non-alias probability ( sometimes called q in other code )
  int    * fAlias; // Alias table           -- could become Index_t ?
};

#endif
