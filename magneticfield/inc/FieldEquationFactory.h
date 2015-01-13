#ifndef FIELDEQUATIONFACTORY_H
#define FIELDEQUATIONFACTORY_H 1

class GUVEquationOfMotion;
class GUVField;

class FieldEquationFactory
{
   static GUVEquationOfMotion*  CreateMagEquation(GUVField *field, int NumEq);
};

#endif
