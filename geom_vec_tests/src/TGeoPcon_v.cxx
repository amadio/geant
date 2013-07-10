#include <iostream>

#include "TGeoPcon_v.h"
#include "TGeoPcon.h"

//_____________________________________________________________________________                                                                                                     
void TGeoPcon_v::Contains_v(const StructOfCoord  & pointi, Bool_t * isin, Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoPcon::Contains(point);
    }
}

//_____________________________________________________________________________

void TGeoPcon_v::Safety_v(const StructOfCoord & pointi, Bool_t in, Double_t * safety, Int_t np ) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      safety[i] = TGeoPcon::Safety(point, in);
    }

}
//_____________________________________________________________________________                                                                                                     
void TGeoPcon_v::DistFromInside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/,const Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoPcon::DistFromInside(point, dir, 3, step[i], 0);
    }
}

//_____________________________________________________________________________                                                                                                    

void TGeoPcon_v::DistFromOutside_v(const StructOfCoord & pointi, const StructOfCoord & diri, Int_t /*iact*/, const  Double_t * step, Double_t * /*safe*/, Double_t * distance , Int_t np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
      Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoPcon::DistFromOutside(point, dir, 3, step[i], 0);
    }

}
