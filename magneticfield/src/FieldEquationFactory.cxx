#include "FieldEquationFactory.h"

// Base types
#include "GUVField.h"
#include "GUVEquationOfMotion.h"

// Concrete Types
#include "TUniformMagField.h"
#include "TMagFieldEquation.h"

GUVEquationOfMotion*
FieldEquationFactory::CreateMagEquation(GUVField *field) // , int NumEq)
{
   const unsigned int NumComp= 6;
   GUVEquationOfMotion *eq= 0;

   TUniformMagField* unifMagFld= dynamic_cast<TUniformMagField*>(field);
   if( unifMagFld != 0 ){
      eq= new TMagFieldEquation<TUniformMagField,NumComp>(unifMagFld); 
   }
   return eq; 
}
