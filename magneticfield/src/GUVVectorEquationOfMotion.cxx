//
// $Id: GUVEquationOfMotion.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// -------------------------------------------------------------------

// #include <iostream>

#include "GUVEquationOfMotion.h"
#include "GUVVectorEquationOfMotion.h"

#include "VcFloatBackend.h"

unsigned int GUVVectorEquationOfMotion::fNumObjectsCreated= 0;
unsigned int GUVVectorEquationOfMotion::fNumObjectsDeleted= 0;

GUVVectorEquationOfMotion::~GUVVectorEquationOfMotion()
{
  CheckDone();
  fNumObjectsDeleted++;
}

void
GUVVectorEquationOfMotion::
EvaluateRhsReturnB( const typename vecgeom::kVc::precision_v  y[],
                          typename vecgeom::kVc::precision_v  dydx[],
                          typename vecgeom::kVc::precision_v  charge,
   vecgeom::Vector3D<typename vecgeom::kVcFloat::precision_v> &Field
                                // vecgeom::Vector3D<Float_v> &Field                    // Tried alternative                    
   ) const
{
   typedef typename vecgeom::kVc::precision_v Double_v;
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
