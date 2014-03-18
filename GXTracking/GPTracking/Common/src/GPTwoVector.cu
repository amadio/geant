#include "GPTwoVector.h"

FQUALIFIER 
GPTwoVector GPTwoVector_create( G4double x, 
				    G4double y)
{
  GPTwoVector v = {x,y};
  return v;
}

FQUALIFIER 
G4double GPTwoVector_x( GPTwoVector v)
{
  return v.x;
}

FQUALIFIER 
G4double GPTwoVector_y( GPTwoVector v)
{
  return v.y;
}

FQUALIFIER 
void GPTwoVector_set(GPTwoVector *v, G4double x, G4double y)
{
  v->x = x;
  v->y = y;
}

FQUALIFIER 
void GPTwoVector_set_x( GPTwoVector *v, G4double x )
{
  v->x = x;
}

FQUALIFIER 
void GPTwoVector_set_y( GPTwoVector *v, G4double y )
{
  v->y = y;
}

FQUALIFIER 
GPTwoVector GPTwoVector_add( GPTwoVector a, GPTwoVector b )
{
  return GPTwoVector_create( a.x+b.x, a.y+b.y);
}

FQUALIFIER 
GPTwoVector GPTwoVector_sub( GPTwoVector a, GPTwoVector b )
{
  return GPTwoVector_create( a.x-b.x, a.y-b.y);
}

FQUALIFIER 
G4double GPTwoVector_mag2( GPTwoVector v )
{
  return v.x*v.x + v.y*v.y;
}

FQUALIFIER 
G4double GPTwoVector_mag( GPTwoVector v )
{
  return sqrt(GPTwoVector_mag2(v));
}

FQUALIFIER 
G4double GPTwoVector_dot( GPTwoVector a, GPTwoVector b )
{
  return a.x*b.x + a.y*b.y;
}

FQUALIFIER 
GPTwoVector GPTwoVector_mult( GPTwoVector a, G4double m )
{
  return GPTwoVector_create( a.x*m, a.y*m);
}

FQUALIFIER 
GPTwoVector GPTwoVector_unit( GPTwoVector v )
{
  G4double mag = GPTwoVector_mag(v);
  if ( mag > 0 ) return GPTwoVector_mult( v, 1.0/mag );
  return v;
}

FQUALIFIER 
G4bool GPTwoVector_equal( GPTwoVector a, GPTwoVector b )
{       
  return a.x == b.x && a.y == b.y ;
}

FQUALIFIER 
G4bool GPTwoVector_nequal( GPTwoVector a, GPTwoVector b )
{       
  return a.x != b.x || a.y != b.y ;
}
