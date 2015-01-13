#include "GUFieldPropagator.h"
#include "TClassicalRK4"

template 
<class T_Field, class T_Stepper, int Neq>
class TFieldPropagator : public GUFieldPropagator
{
   TFieldPropagator();
   ~TFieldPropagator(){ delete fEquation; delete fStepper; delete fDriver;}

   GUVEquationOfMotion*   fEquation;
   GUVIntegrationStepper* fStepper;
   GUIntegrationDriver*   fDriver;
};

template 
<class T_Field, class T_Stepper, Neq>
TFieldPropagator::TFieldPropagator(typename T_Field* field)
{
  // Must create the Driver, Stepper and Equation ??
  TMagFieldEquation*  pEquation = new TMagFieldEquation<field, Neq>;
  fEquation = pEquation;
  // GUVIntegrationStepper* 
  fStepper = new TClassicalRK4<pEquation,Neq>;
  fDriver  = new TIntegrationDriver<fStepper>;
}