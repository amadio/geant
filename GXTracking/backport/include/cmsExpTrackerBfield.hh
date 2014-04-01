#ifndef cmsExpTrackerBfield_H
#define cmsExpTrackerBfield_H 1

// B-field in the CMS Tracker volume - based on the TOSCA computation
//   In :   x[3]: coordinates (m)
//   Out:   B[3]: Bx,By,Bz    (T) (getBxyz)
//   Out:   B[3]: Br,Bf,Bz    (T) (getBrfz)
//
// Valid for r<1.15m and |z|<2.80m, originally written by V. Karimaki
// Modified for the standalone application with CMS detector geometry

//#include "globals.hh"
#include <string>

class cmsExpTrackerBfield
{
  public:

  cmsExpTrackerBfield ();  
  ~cmsExpTrackerBfield() {}

  // B out in cartesian
  void getBxyz(G4double const *x, G4double *Bxyz) const; 
  // B out in cylindrical
  void getBrfz(G4double const *x, G4double *Brfz) const;

  private:
  G4double prm[9];
  G4double ap2, hb0, hlova, ainv, coeff;
  void Bcyl(G4double r, G4double z, G4double *Bw) const;
};

#endif
