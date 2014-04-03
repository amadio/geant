#ifndef GPThreeVector_h
#define GPThreeVector_h 1

#include "GPTypeDef.h"
#include "GPGeomdefs.h"

struct GPThreeVector
{
  double x;
  double y;
  double z;
};

extern "C" {

FQUALIFIER
GPThreeVector GPThreeVector_create( double x, 
				    double y, 
				    double z );
FQUALIFIER 
G4double GPThreeVector_x( GPThreeVector v);

  FQUALIFIER 
G4double GPThreeVector_y( GPThreeVector v);

FQUALIFIER 
G4double GPThreeVector_z( GPThreeVector v);

FQUALIFIER 
void GPThreeVector_set(GPThreeVector *v, G4double x, G4double y, G4double z );

FQUALIFIER 
void GPThreeVector_set_x( GPThreeVector *v, G4double x );

FQUALIFIER 
void GPThreeVector_set_y( GPThreeVector *v, G4double y );

FQUALIFIER 
GPThreeVector GPThreeVector_add( GPThreeVector a, GPThreeVector b );

FQUALIFIER 
GPThreeVector GPThreeVector_sub( GPThreeVector a, GPThreeVector b );

FQUALIFIER 
double GPThreeVector_mag2( GPThreeVector v );

FQUALIFIER 
double GPThreeVector_mag( GPThreeVector v );

FQUALIFIER 
double GPThreeVector_dot( GPThreeVector a, GPThreeVector b );

FQUALIFIER 
GPThreeVector GPThreeVector_cross( GPThreeVector a, GPThreeVector p );

FQUALIFIER 
GPThreeVector GPThreeVector_mult( GPThreeVector a, double m );

FQUALIFIER 
GPThreeVector GPThreeVector_unit( GPThreeVector v );

FQUALIFIER 
G4bool GPThreeVector_equal( GPThreeVector a, GPThreeVector b );

FQUALIFIER 
G4bool GPThreeVector_nequal( GPThreeVector a, GPThreeVector b );

FQUALIFIER
void GPThreeVector_rotateUz( GPThreeVector* This, GPThreeVector reference );

FQUALIFIER
GPThreeVector GPThreeVector_saxpy(G4double a, GPThreeVector x, GPThreeVector y);

FQUALIFIER
G4double GPThreeVector_coord( GPThreeVector v, EAxis axis );

FQUALIFIER
void GPThreeVector_set_coord( GPThreeVector *v, EAxis axis, G4double val );

}

#endif
