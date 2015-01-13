//
// $Id: GUVEquationOfMotion.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// -------------------------------------------------------------------

#include "GUVEquationOfMotion.h"

GUVEquationOfMotion::~GUVEquationOfMotion()
{}

void 
GUVEquationOfMotion::
EvaluateRhsReturnB( const double y[],
     				 double dydx[],
	         			 double charge, 
                          double Field[] 
                  ) const
{
   double  PositionAndTime[4];
   // Position
   PositionAndTime[0] = y[0];
   PositionAndTime[1] = y[1];
   PositionAndTime[2] = y[2];
   // Global Time
   PositionAndTime[3] = y[7];  // See GUVFieldTrack::LoadFromArray

   GetFieldValue(PositionAndTime, Field) ;
   EvaluateRhsGivenB( y, Field, charge, dydx );
}

// #if  HELP_THE_COMPILER
// void 
// GUVEquationOfMotion::doNothing()
// {
// }
// #endif
