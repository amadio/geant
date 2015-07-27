#ifndef ROOT_TGeoPcon_v
#define ROOT_TGeoPcon_v 

#ifndef ROOT_TGeoShape
#include "TGeoShape.h"
#endif
#include "TGeoPcon.h"

#ifdef VEC_EXTENSIONS
#include "Vc/vector.h"
#endif


///////////////////////////////////////////////////////////////////////////////

#include "PointStruct.h" // for SOA data

class TGeoPcon_v: public TGeoPcon
{
public:
  virtual void Contains_v(const StructOfCoord & pointi, Bool_t * isin, Int_t np) const; 
  virtual void Safety_v(const StructOfCoord &pointi, Bool_t in, double *safety, Int_t np ) const;         
  virtual void DistFromInside_v(StructOfCoord const  & pointi, const StructOfCoord & diri, Int_t iact, const double * step, double *safe, double * distance , Int_t np) const; 
  virtual void DistFromOutside_v(const StructOfCoord & pointi, const  StructOfCoord & diri, Int_t iact, const double * step, double *safe, double * distance , Int_t np) const; 

};
#endif
