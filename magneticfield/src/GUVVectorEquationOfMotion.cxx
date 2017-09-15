//
// $Id: GUVEquationOfMotion.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// -------------------------------------------------------------------

// #include <iostream>

#include "GUVEquationOfMotion.h"
#include "GUVVectorEquationOfMotion.h"

unsigned int GUVVectorEquationOfMotion::fNumObjectsCreated= 0;
unsigned int GUVVectorEquationOfMotion::fNumObjectsDeleted= 0;

GUVVectorEquationOfMotion::~GUVVectorEquationOfMotion()
{
  fNumObjectsDeleted++;
}

void
GUVVectorEquationOfMotion::
EvaluateRhsReturnB( const Double_v  y[],
                          Double_v  dydx[],
                          Double_v  charge,
                          vecgeom::Vector3D<Double_v> &Field
   ) const
{
   Double_v  PositionAndTime[4];
   PositionAndTime[0] = y[0];
   PositionAndTime[1] = y[1];
   PositionAndTime[2] = y[2];
   PositionAndTime[3] = y[7];  

   GetFieldValue( PositionAndTime, Field) ;
   EvaluateRhsGivenB( y, Field, charge, dydx );
}

std::ostream&  operator<<( std::ostream& os, const GUVVectorEquationOfMotion& eq)
{
   os << " Equation of Motion # " << eq.GetId()
      << "   field ptr= "  << eq.GetFieldObj() << "  Initialised= " << eq.Initialised()
      << std::endl;
   os << "  Total # of E-of-M = " << GUVVectorEquationOfMotion::GetNumCreated()
      << " live= " << GUVVectorEquationOfMotion::GetNumLive() << std::endl;
   return os;
}
