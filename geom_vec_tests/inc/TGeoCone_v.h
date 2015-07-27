#ifndef ROOT_TGeoPcon_v
#define ROOT_TGeoPcon_v

#ifndef ROOT_TGeoShape
#include "TGeoShape.h"
#endif
#include "TGeoCone.h"

#ifdef VEC_EXTENSIONS
#include "Vc/vector.h"
#endif



///////////////////////////////////////////////////////////////////////////////

#include "PointStruct.h" // for SOA data

class TGeoCone_v: public TGeoCone /*public TGeoConeSeg*/
{
public:
    virtual void Contains_v(const StructOfCoord & pointi, Bool_t * isin, Int_t np) const;
    
    virtual void Safety_v(const StructOfCoord &pointi, Bool_t in, double *safety, Int_t np ) const;
    virtual void DistFromInside_v(StructOfCoord const  & pointi, const StructOfCoord & diri, Int_t iact, const double * step, double *safe, double * distance , Int_t np) const;
    virtual void DistFromOutside_v(const StructOfCoord & pointi, const  StructOfCoord & diri, Int_t iact, const double * step, double *safe, double * distance , Int_t np) const;
    //virtual void SafetySeg_v(Vc::double_v r, Vc::double_v z, double r1,double z1, double r2, double z2, Bool_t outer, double * safetySeg, Int_t np ) const;
    
    
};

class TGeoConeSeg_v: public TGeoConeSeg  /*public TGeoCone_v*/
{
public:
    virtual void Contains_v(const StructOfCoord & pointi, Bool_t * isin, Int_t np) const;
#ifdef VEC_EXTENSIONS
private:
    virtual void Contains_v4( Vc::double_v const & x , Vc::double_v const & y, Vc::double_v const & z,  Vc::double_m &c1);
public:
#endif
    virtual void Safety_v(const StructOfCoord &pointi, Bool_t in, double *safety, Int_t np ) const;
    virtual void DistFromInside_v(StructOfCoord const  & pointi, const StructOfCoord & diri, Int_t iact, const double * step, double *safe, double * distance , Int_t np) const;
    virtual void DistFromOutside_v(const StructOfCoord & pointi, const  StructOfCoord & diri, Int_t iact, const double * step, double *safe, double * distance , Int_t np) const;
    
    
};

#endif

