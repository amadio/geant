#ifndef GPTwoVector_H
#define GPTwoVector_H 1

#include "GPTypeDef.h"

struct GPTwoVector
{
  G4double x;
  G4double y;
};

extern "C" {

FQUALIFIER 
GPTwoVector GPTwoVector_create( G4double x, 
				G4double y);

FQUALIFIER 
G4double GPTwoVector_x( GPTwoVector v);

FQUALIFIER 
G4double GPTwoVector_y( GPTwoVector v);

FQUALIFIER 
void GPTwoVector_set(GPTwoVector *v, G4double x, G4double y);

FQUALIFIER 
void GPTwoVector_set_x( GPTwoVector *v, G4double x );

FQUALIFIER 
void GPTwoVector_set_y( GPTwoVector *v, G4double y );

FQUALIFIER 
GPTwoVector GPTwoVector_add( GPTwoVector a, GPTwoVector b );

FQUALIFIER 
GPTwoVector GPTwoVector_sub( GPTwoVector a, GPTwoVector b );

FQUALIFIER 
G4double GPTwoVector_mag2( GPTwoVector v );

FQUALIFIER 
G4double GPTwoVector_mag( GPTwoVector v );

FQUALIFIER 
G4double GPTwoVector_dot( GPTwoVector a, GPTwoVector b );

FQUALIFIER 
GPTwoVector GPTwoVector_mult( GPTwoVector a, G4double m );

FQUALIFIER 
GPTwoVector GPTwoVector_unit( GPTwoVector v );

FQUALIFIER 
G4bool GPTwoVector_equal( GPTwoVector a, GPTwoVector b );

FQUALIFIER 
G4bool GPTwoVector_nequal( GPTwoVector a, GPTwoVector b );

}

#endif
