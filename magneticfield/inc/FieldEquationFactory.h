#ifndef FIELDEQUATIONFACTORY_H
#define FIELDEQUATIONFACTORY_H 1

// Base types
// #include "GUVField.h"
#include "GUVEquationOfMotion.h"

// Concrete Types
// #include "TUniformMagField.h"
#include "TMagFieldEquation.h"

class FieldEquationFactory
{
   public:
     static const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

     template<typename FieldType>
        static
        TMagFieldEquation<FieldType,Nposmom> *   // GUVEquationOfMotion*
               CreateMagEquation(FieldType *field);
       //  Create an equation given a field type

     template<typename FieldType>
        static
        TMagFieldEquation<FieldType,Nposmom> *       
               CreateMagEquation(FieldType &field);
      //  Similar for a field reference 
};

template<typename FieldType>
TMagFieldEquation<FieldType,FieldEquationFactory::Nposmom> *   // GUVEquationOfMotion*
FieldEquationFactory::CreateMagEquation(FieldType *pField)
{
    return new TMagFieldEquation<FieldType, Nposmom>(pField);
}

template<typename FieldType>
TMagFieldEquation<FieldType,FieldEquationFactory::Nposmom> *   // GUVEquationOfMotion*
FieldEquationFactory::CreateMagEquation(FieldType &field)
{
    return new TMagFieldEquation<FieldType, Nposmom>(&field);
}
#endif
