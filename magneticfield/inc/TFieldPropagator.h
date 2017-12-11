#include "GUFieldPropagator.h"
#include "TClassicalRK4.h"

template 
<class T_Equation, class T_Stepper, int Neq>
class TFieldPropagator : public GUFieldPropagator
{
   TFieldPropagator(T_Equation equation, T_Stepper stepper, driver);
   ~TFieldPropagator(){ delete fEquation; delete fStepper; delete fDriver;}

   VScalarEquationOfMotion*   fEquation;
   VScalarIntegrationStepper* fStepper;
   GUIntegrationDriver*   fDriver;
};

template 
<class T_Field, class T_Stepper, int Neq>
TFieldPropagator::TFieldPropagator(typename T_Field* field)
{
  // Must create the Driver, Stepper and Equation ??
  TMagFieldEquation*  pEquation = new TMagFieldEquation<field, Neq>;
  fEquation = pEquation;
  // VScalarIntegrationStepper* 
  fStepper = new TClassicalRK4<pEquation,Neq>;
  fDriver  = new TIntegrationDriver<fStepper>;
}
