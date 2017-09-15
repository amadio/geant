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
   fNumObjectsDeleted++;
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
