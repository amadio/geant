#ifndef GXSimpleEcal_HH
#define GXSimpleEcal_HH

#include "GXUserGeometry.h"

class GXSimpleEcal : public GXUserGeometry
{
public:
  GXSimpleEcal( int nphi=3, int nz=4, double density=8.28*g/cm3);
   
  inline double getScale() const { return 1.0; }

private:
   void create_impl();
   
private:
   int fNphi;
   int fNz;
   double fDensity;
   
};

#endif
