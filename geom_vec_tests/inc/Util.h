#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include "TRandom.h"

class Util{
  public:
  static void samplePoint( double *point, double dx, double dy, double dz, double scale )
  {
    point[0]=scale*(1-2.*gRandom->Rndm())*dx;
    point[1]=scale*(1-2.*gRandom->Rndm())*dy;
    point[2]=scale*(1-2.*gRandom->Rndm())*dz;
  }
 
  static void samplePoint( double *point, double dx, double dy, double dz, double const * origin, double scale )
  {
    point[0]=origin[0]+scale*(1-2.*gRandom->Rndm())*dx;
    point[1]=origin[1]+scale*(1-2.*gRandom->Rndm())*dy;
    point[2]=origin[2]+scale*(1-2.*gRandom->Rndm())*dz;
  }
  

  // creates random normalized vectors 
  static void sampleDir( double * dir )
  {
    dir[0]=(1-2.*gRandom->Rndm());
    dir[1]=(1-2.*gRandom->Rndm());
    dir[2]=(1-2.*gRandom->Rndm());
    double inversenorm=1./sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
    dir[0]*=inversenorm;
    dir[1]*=inversenorm;
    dir[2]*=inversenorm;
  }


};


#endif// UTIL_H
