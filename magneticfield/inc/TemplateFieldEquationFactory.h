#ifndef TEMPLATEFIELDEQUATIONFACTORY_H
#define TEMPLATEFIELDEQUATIONFACTORY_H 1


#include "TemplateGUVEquationOfMotion.h"

#include "TemplateTMagFieldEquation.h"

template <class Backend>
class TemplateFieldEquationFactory
{
  public:
    static const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

    template<typename FieldType>
    static
    TemplateTMagFieldEquation<Backend,FieldType,Nposmom> * CreateMagEquation(FieldType *field);
       //  Create an equation given a field type
};

template <class Backend>
template<typename FieldType>
TemplateTMagFieldEquation<Backend, FieldType,TemplateFieldEquationFactory<Backend>::Nposmom> *   // GUVEquationOfMotion*
TemplateFieldEquationFactory<Backend>
  ::CreateMagEquation(FieldType *field)
{
  return new TemplateTMagFieldEquation<Backend, FieldType, Nposmom>(field);
}
#endif
