typedef double  G4double;
typedef int G4int;

#include "cmsExpTrackerBfield.hh"
#include "cmsExpNameSpace.hh"
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "G4SystemOfUnits.hh"

namespace cmsExp {

//beginning of the cmsExp name space

void ffunkti(G4double u, G4double *ff)
{
  // Function and its 3 derivatives
  G4double a = 1.0/(1.0+u*u);
  G4double b = std::sqrt(a);

  ff[0]=u*b;
  ff[1]=a*b;
  ff[2]=-3.0*u*a*ff[1];
  ff[3]=a*ff[2]*((1.0/u)-4.0*u);
}

//end of the cmxExp name space
}

cmsExpTrackerBfield::cmsExpTrackerBfield() {

  for (G4int i=0; i<9; i++) prm[i] = cmsExp::Bp_3_8T[i];

  if (!prm[0]) {
    std::cerr << "cmsExpTrackerBfield::cmsExpTrackerBfield(), BadParameters: Defined keys are: 3_8T";
    exit (1);
  }

  ap2=4*prm[0]*prm[0]/(prm[1]*prm[1]);  
  hb0=0.5*prm[2]*std::sqrt(1.0+ap2);
  hlova=1/std::sqrt(ap2);
  ainv=2*hlova/prm[1];
  coeff=1/(prm[8]*prm[8]);

}

void  cmsExpTrackerBfield::Bcyl(G4double r, G4double z, G4double *Bw) const
{
  // convert to m (cms magnetic field parameterization)
  r *=0.001;
  z *=0.001;

  z-=prm[3];                    // max Bz point is shifted in z

  G4double az=std::abs(z);
  G4double zainv=z*ainv;
  G4double u=hlova-zainv;
  G4double v=hlova+zainv;
  G4double fu[4],gv[4];

  cmsExp::ffunkti(u,fu);
  cmsExp::ffunkti(v,gv);

  G4double rat=0.5*r*ainv;
  G4double rat2=rat*rat;
  Bw[0]=hb0*rat*(fu[1]-gv[1]-(fu[3]-gv[3])*rat2*0.5);
  Bw[1]=0;
  Bw[2]=hb0*(fu[0]+gv[0]-(fu[2]+gv[2])*rat2);
  G4double corBr= prm[4]*r*z*(az-prm[5])*(az-prm[5]);
  G4double corBz=-prm[6]*(exp(-(z-prm[7])*(z-prm[7])*coeff) +
			  exp(-(z+prm[7])*(z+prm[7])*coeff) );
  Bw[0]+=corBr;
  Bw[2]+=corBz;
}

void cmsExpTrackerBfield::getBxyz(G4double const *x, G4double *Bxyz)  const 
{
   G4double Bw[3];

  G4double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  Bcyl(r, x[2], Bw);
  G4double rinv=(r>0) ? 1/r:0;
  Bxyz[0]=tesla*Bw[0]*x[0]*rinv;
  Bxyz[1]=tesla*Bw[0]*x[1]*rinv;
  Bxyz[2]=tesla*Bw[2];

}
