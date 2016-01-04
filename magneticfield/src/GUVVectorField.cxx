//
// First implementation class for GUVVectorField 
// -------------------------------------------------------------------

#include "GUVField.h"
#include "GUVVectorField.h"
#include <iostream>

 
GUVVectorField::~GUVVectorField()
{
}

//Confirm about commenting this function. Assumed same functionality derived from GUVField
GUVVectorField& GUVVectorField::operator = (const GUVVectorField &p)
{
  
   if (&p != this){

    //line below if inheriting from GUVField. Comment 2nd line in that case
    // this->GUVField::operator=(p);
   
     this->fChangesEnergy= p.fChangesEnergy;   
   }
   return *this;
}


GUVVectorField* GUVVectorField::Clone() const
{
    std::cout << "Derived class does not implement cloning,\n"
              << "but Clone method called.\n"
              << "Cannot continue;" << std::endl;
    exit(1); 
    return NULL;
}
// ------------------------------------------------------------------------
