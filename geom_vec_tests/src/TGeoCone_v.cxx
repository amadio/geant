#include <iostream>
#include "TGeoCone_v.h"


#ifndef VEC_EXTENSIONS
//_____________________________________________________________________________                                                                                                     
void TGeoCone_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoCone::Contains(point);
    }
}
#else
// PUT VC CODE OR THE LIKE HERE
//_____________________________________________________________________________                                                                                                     
void TGeoCone_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoCone::Contains(point);
    }
}
#endif


//_____________________________________________________________________________

void TGeoCone_v::Safety_v(const StructOfCoord & pointi, Bool_t in, Double_t * safety, Int_t np ) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      safety[i] = TGeoCone::Safety(point, in);
    }
}

//_____________________________________________________________________________                                                                                                     
void TGeoCone_v::DistFromInside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/,const Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoCone::DistFromInside(point, dir, 3, step[i], 0);
    }
}

//_____________________________________________________________________________                                                                                                    

void TGeoCone_v::DistFromOutside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/, const  Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoCone::DistFromOutside(point, dir, 3, step[i], 0);
    }

}


//_____________________________________________________________________________                                                                                                     
void TGeoConeSeg_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoConeSeg::Contains(point);
    }
}

//_____________________________________________________________________________

void TGeoConeSeg_v::Safety_v(const StructOfCoord & pointi, Bool_t in, Double_t * safety, Int_t np ) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      safety[i] = TGeoConeSeg::Safety(point, in);
    }

}
//_____________________________________________________________________________                                                                                                     
void TGeoConeSeg_v::DistFromInside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/,const Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoConeSeg::DistFromInside(point, dir, 3, step[i], 0);
    }
}

//_____________________________________________________________________________                                                                                                    

void TGeoConeSeg_v::DistFromOutside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/, const  Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoConeSeg::DistFromOutside(point, dir, 3, step[i], 0);
    }

}
