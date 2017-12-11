#ifndef FIELDEQUATIONFACTORY_H
#define FIELDEQUATIONFACTORY_H 1

// Base types
// #include "GUVField.h"
#include "VScalarEquationOfMotion.h"

// Concrete Types
// #include "ScalarUniformMagField.h"
#include "ScalarFieldEquation.h"

class FieldEquationFactory
{
   public:
     static const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

     template<typename FieldType>
        static
        ScalarFieldEquation<FieldType,Nposmom> *   // VScalarEquationOfMotion*
               CreateMagEquation(FieldType *field);
       //  Create an equation given a field type

     template<typename FieldType>
        static
        ScalarFieldEquation<FieldType,Nposmom> *       
               CreateMagEquation(FieldType &field);
      //  Similar for a field reference 
};

template<typename FieldType>
ScalarFieldEquation<FieldType,FieldEquationFactory::Nposmom> *   // VScalarEquationOfMotion*
FieldEquationFactory::CreateMagEquation(FieldType *pField)
{
    return new ScalarFieldEquation<FieldType, Nposmom>(pField);
}

template<typename FieldType>
ScalarFieldEquation<FieldType,FieldEquationFactory::Nposmom> *   // VScalarEquationOfMotion*
FieldEquationFactory::CreateMagEquation(FieldType &field)
{
    return new ScalarFieldEquation<FieldType, Nposmom>(&field);
}
#endif
