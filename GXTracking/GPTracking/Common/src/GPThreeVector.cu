#include <math.h>
#include "GPThreeVector.h"

FQUALIFIER 
GPThreeVector GPThreeVector_create( double x, 
				    double y, 
				    double z )
{
  GPThreeVector v = {x,y,z};
  return v;
}

FQUALIFIER 
G4double GPThreeVector_x( GPThreeVector v)
{
  return v.x;
}

FQUALIFIER 
G4double GPThreeVector_y( GPThreeVector v)
{
  return v.y;
}

FQUALIFIER 
G4double GPThreeVector_z( GPThreeVector v)
{
  return v.z;
}

FQUALIFIER 
void GPThreeVector_set(GPThreeVector *v, G4double x, G4double y, G4double z )
{
  v->x = x;
  v->y = y;
  v->z = z;
}

FQUALIFIER 
void GPThreeVector_set_x( GPThreeVector *v, G4double x )
{
  v->x = x;
}

FQUALIFIER 
void GPThreeVector_set_y( GPThreeVector *v, G4double y )
{
  v->y = y;
}

FQUALIFIER 
void GPThreeVector_set_z( GPThreeVector *v, G4double z )
{
  v->z = z;
}

FQUALIFIER 
GPThreeVector GPThreeVector_add( GPThreeVector a, GPThreeVector b )
{
  return GPThreeVector_create( a.x+b.x, a.y+b.y, a.z+b.z );
}

FQUALIFIER 
GPThreeVector GPThreeVector_sub( GPThreeVector a, GPThreeVector b )
{
  return GPThreeVector_create( a.x-b.x, a.y-b.y, a.z-b.z );
}

FQUALIFIER 
double GPThreeVector_mag2( GPThreeVector v )
{
  return v.x*v.x + v.y*v.y + v.z*v.z;
}

FQUALIFIER 
double GPThreeVector_mag( GPThreeVector v )
{
  return sqrt(GPThreeVector_mag2(v));
}

FQUALIFIER 
double GPThreeVector_dot( GPThreeVector a, GPThreeVector b )
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

FQUALIFIER 
GPThreeVector GPThreeVector_cross( GPThreeVector a, GPThreeVector p )
{
  return GPThreeVector_create( a.y*p.z-p.y*a.z,
			       a.z*p.x-p.z*a.x,
			       a.x*p.y-p.x*a.y );
}

FQUALIFIER 
GPThreeVector GPThreeVector_mult( GPThreeVector a, double m )
{
  return GPThreeVector_create( a.x*m, a.y*m, a.z*m );
}

FQUALIFIER 
GPThreeVector GPThreeVector_unit( GPThreeVector v )
{
  double mag = GPThreeVector_mag(v);
  if ( mag > 0 ) return GPThreeVector_mult( v, 1.0/mag );
  return v;
}

FQUALIFIER 
G4bool GPThreeVector_equal( GPThreeVector a, GPThreeVector b )
{       
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

FQUALIFIER 
G4bool GPThreeVector_nequal( GPThreeVector a, GPThreeVector b )
{       
  return a.x != b.x || a.y != b.y || a.z != b.z;
}

FQUALIFIER
void GPThreeVector_rotateUz( GPThreeVector* This, GPThreeVector reference )
{
  // NewUzVector must be normalized !
  GPThreeVector NewUzVector = GPThreeVector_unit(reference);

  double u1 = NewUzVector.x;
  double u2 = NewUzVector.y;
  double u3 = NewUzVector.z;
  double up = u1*u1 + u2*u2;

  if (up>0) {
      up = sqrt(up);
      double px = This->x,  py = This->y,  pz = This->z;
      This->x = (u1*u3*px - u2*py)/up + u1*pz;
      This->y = (u2*u3*px + u1*py)/up + u2*pz;
      This->z =    -up*px +             u3*pz;
    }
  else if (u3 < 0.) { This->x = -This->x; This->z = -This->z; } //phi=0 teta=pi
  else {
    ;
  }

}

FQUALIFIER
GPThreeVector GPThreeVector_saxpy(G4double a, GPThreeVector x, GPThreeVector y)
{
  return GPThreeVector_create(a*x.x + y.x,
			      a*x.y + y.y,
			      a*x.z + y.z );
}

FQUALIFIER
G4double GPThreeVector_coord( GPThreeVector v, EAxis axis )
{
        switch( axis )
        {
        case kXAxis: return v.x;
        case kYAxis: return v.y;
        case kZAxis: return v.z;
        default:
	  //                myAssert(false);
                return 0;
        }
}

FQUALIFIER
void GPThreeVector_set_coord( GPThreeVector *v, EAxis axis, G4double val )
{
        switch( axis )
        {
        case kXAxis: v->x = val; break;
        case kYAxis: v->y = val; break;
        case kZAxis: v->z = val; break;
        default:
	  //                myAssert(false);
                break;
        }
}
