#ifndef FIELDEQUATIONFACTORY_H
#define FIELDEQUATIONFACTORY_H 1

// Base types
// #include "VVectorField.h"
#include "VScalarEquationOfMotion.h"

// Concrete Types
// #include "ScalarUniformMagField.h"
#include "ScalarMagFieldEquation.h"

class FieldEquationFactory
{
   public:
     static const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

     template<typename FieldType>
        static
        ScalarMagFieldEquation<FieldType,Nposmom> *   // VScalarEquationOfMotion*
               CreateMagEquation(FieldType *field);
       //  Create an equation given a field type

     template<typename FieldType>
        static
        ScalarMagFieldEquation<FieldType,Nposmom> *       
               CreateMagEquation(FieldType &field);
      //  Similar for a field reference 
};

template<typename FieldType>
ScalarMagFieldEquation<FieldType,FieldEquationFactory::Nposmom> *   // VScalarEquationOfMotion*
FieldEquationFactory::CreateMagEquation(FieldType *pField)
{
    return new ScalarMagFieldEquation<FieldType, Nposmom>(pField);
}

template<typename FieldType>
ScalarMagFieldEquation<FieldType,FieldEquationFactory::Nposmom> *   // VScalarEquationOfMotion*
FieldEquationFactory::CreateMagEquation(FieldType &field)
{
    return new ScalarMagFieldEquation<FieldType, Nposmom>(&field);
}
#endif
