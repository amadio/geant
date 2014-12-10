//
// First implementation class for GUVField 
// -------------------------------------------------------------------

#include "GUVField.h"
#include <iostream>

GUVField::GUVField( )
{
}
 
GUVField::~GUVField()
{
}

GUVField& GUVField::operator = (const GUVField &p)
{
   if (&p == this) return *this;
   return *this;
}

GUVField::GUVField (const GUVField &p)
{
}

GUVField* GUVField::Clone() const
{
   // GUExceptionDescription msg;
    std::cout << "Derived class does not implement cloning,\n"
              << "but Clone method called.\n"
              << "Cannot continue;" << std::endl;
    exit(1); 
    // GUException("GUVField::Clone", "GeomField004", FatalException, msg );
    return NULL;
}
// ------------------------------------------------------------------------
