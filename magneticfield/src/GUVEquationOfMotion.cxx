//
// $Id: GUVEquationOfMotion.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// -------------------------------------------------------------------

// #include <iostream>

#include "GUVEquationOfMotion.h"

unsigned int GUVEquationOfMotion::fNumObjectsCreated= 0;
unsigned int GUVEquationOfMotion::fNumObjectsDeleted= 0;

GUVEquationOfMotion::~GUVEquationOfMotion()
{
   CheckDone();
   // std::cout << " Destructing Equation " << this << " info= " << *this << std::endl;   
   fNumObjectsDeleted++;
   // To help ensure that clients call InformDone() - ie. clear
}

void
GUVEquationOfMotion::
EvaluateRhsReturnB( const double           y[],
                          double          dydx[],
                       // double          charge,
                 vecgeom::Vector3D<float> &Field
                  ) const
{
   double  PositionAndTime[4];
   // Position
   PositionAndTime[0] = y[0];
   PositionAndTime[1] = y[1];
   PositionAndTime[2] = y[2];
   // Global Time
   PositionAndTime[3] = y[7];  // See GUVFieldTrack::LoadFromArray

   GetFieldValue( PositionAndTime, Field) ;
   EvaluateRhsGivenB( y, Field, /*charge,*/ dydx );
}

std::ostream&  operator<<( std::ostream& os, const GUVEquationOfMotion& eq)
{
   os << " Equation of Motion # " << eq.GetId()
      << "   field ptr= "  << eq.GetFieldObj() << "  Initialised= " << eq.Initialised()
      << std::endl;
   os << "  Total # of E-of-M = " << GUVEquationOfMotion::GetNumCreated()
      << " live= " << GUVEquationOfMotion::GetNumLive() << std::endl;
   return os;
}
