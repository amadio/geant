#include "FieldEquationFactory.h"

#include "GUVEquationOfMotion.h"
#include "TMagFieldEquation.h"
#include "GUVField.h"

GUVEquationOfMotion*
FieldEquationFactory::CreateMagEquation(GUVField *field, int NumEq)
{
   GUVEquationOfMotion *eq= 0;
   if( dynamic_cast<GUUniformMagField>(field) != 0 ) {
      eq= TMagFieldEquation<GUUniformMagField,NumEq>(field); 
   }
   return eq; 
}
