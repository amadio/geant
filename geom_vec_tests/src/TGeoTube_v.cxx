#include <iostream>
#include "TGeoTube_v.h"


#ifndef VEC_EXTENSIONS
//_____________________________________________________________________________                                                                                                     
void TGeoTube_v::Contains_v(const StructOfCoord  & pointi, bool * isin, int np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoTube::Contains(point);
    }
}
#else
// PUT VC CODE OR THE LIKE HERE
void TGeoTube_v::Contains_v(const StructOfCoord  & pointi, bool * isin, int np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoTube::Contains(point);
    }
}
#endif


//_____________________________________________________________________________

void TGeoTube_v::Safety_v(const StructOfCoord & pointi, bool in, double * safety, int np ) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      safety[i] = TGeoTube::Safety(point, in);
    }
}

//_____________________________________________________________________________                                                                                                     
void TGeoTube_v::DistFromInside_v(const StructOfCoord & pointi, const StructOfCoord & diri, int /*iact*/,const double * step, double * /*safe*/, double * distance , int np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      double dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoTube::DistFromInside(point, dir, 3, step[i], 0);
    }
}

//_____________________________________________________________________________                                                                                                    

void TGeoTube_v::DistFromOutside_v(const StructOfCoord & pointi, const StructOfCoord & diri, int /*iact*/, const  double * step, double * /*safe*/, double * distance , int np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
      double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      double dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoTube::DistFromOutside(point, dir, 3, step[i], 0);
    }

}


//_____________________________________________________________________________                                                                                                     
void TGeoTubeSeg_v::Contains_v(const StructOfCoord  & pointi, bool * isin, int np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      isin[i] = TGeoTubeSeg::Contains(point);
    }
}

//_____________________________________________________________________________

void TGeoTubeSeg_v::Safety_v(const StructOfCoord & pointi, bool in, double * safety, int np ) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      safety[i] = TGeoTubeSeg::Safety(point, in);
    }

}
//_____________________________________________________________________________                                                                                                     
void TGeoTubeSeg_v::DistFromInside_v(const StructOfCoord & pointi, const StructOfCoord & diri, int /*iact*/,const double * step, double * /*safe*/, double * distance , int np) const
{
  for(unsigned int i = 0; i < np; i++)
    {
      double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      double dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoTubeSeg::DistFromInside(point, dir, 3, step[i], 0);
    }
}

//_____________________________________________________________________________                                                                                                    

void TGeoTubeSeg_v::DistFromOutside_v(const StructOfCoord & pointi, const StructOfCoord & diri, int /*iact*/, const  double * step, double * /*safe*/, double * distance , int np) const
{
    for(unsigned int i = 0; i < np; i++)
    {
      double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
      double dir[3] = {diri.x[i], diri.y[i], diri.z[i]};
      
      distance[i] = TGeoTubeSeg::DistFromOutside(point, dir, 3, step[i], 0);
    }

}
