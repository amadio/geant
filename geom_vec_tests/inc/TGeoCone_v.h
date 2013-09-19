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
    
    virtual void Safety_v(const StructOfCoord &pointi, Bool_t in, Double_t *safety, Int_t np ) const;
    virtual void DistFromInside_v(StructOfCoord const  & pointi, const StructOfCoord & diri, Int_t iact, const Double_t * step, Double_t *safe, Double_t * distance , Int_t np) const;
    virtual void DistFromOutside_v(const StructOfCoord & pointi, const  StructOfCoord & diri, Int_t iact, const Double_t * step, Double_t *safe, Double_t * distance , Int_t np) const;
    //virtual void SafetySeg_v(Vc::double_v r, Vc::double_v z, Double_t r1,Double_t z1, Double_t r2, Double_t z2, Bool_t outer, Double_t * safetySeg, Int_t np ) const;
    
    
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
    virtual void Safety_v(const StructOfCoord &pointi, Bool_t in, Double_t *safety, Int_t np ) const;
    virtual void DistFromInside_v(StructOfCoord const  & pointi, const StructOfCoord & diri, Int_t iact, const Double_t * step, Double_t *safe, Double_t * distance , Int_t np) const;
    virtual void DistFromOutside_v(const StructOfCoord & pointi, const  StructOfCoord & diri, Int_t iact, const Double_t * step, Double_t *safe, Double_t * distance , Int_t np) const;
    
    
};

#endif

