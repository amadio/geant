//
// $Id: VScalarEquationOfMotion.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// -------------------------------------------------------------------

// #include <iostream>

#include "VScalarEquationOfMotion.h"

unsigned int VScalarEquationOfMotion::fNumObjectsCreated= 0;
unsigned int VScalarEquationOfMotion::fNumObjectsDeleted= 0;

VScalarEquationOfMotion::~VScalarEquationOfMotion()
{
   fNumObjectsDeleted++;
}


std::ostream&  operator<<( std::ostream& os, const VScalarEquationOfMotion& eq)
{
   os << " Equation of Motion # " << eq.GetId()
      << "   field ptr= "  << eq.GetFieldObj() << "  Initialised= " << eq.Initialised()
      << std::endl;
   os << "  Total # of E-of-M = " << VScalarEquationOfMotion::GetNumCreated()
      << " live= " << VScalarEquationOfMotion::GetNumLive() << std::endl;
   return os;
}
