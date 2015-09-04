//
//
// Derived from G4MagIntegrationStepper class 
//
// --------------------------------------------------------------------

#include "GUVIntegrationStepper.h"

// Constructor for stepper abstract base class. 
// 

GUVIntegrationStepper::GUVIntegrationStepper(GUVEquationOfMotion* equation,
					     unsigned int num_integration_vars,
                                             unsigned int integrationOrder,
					     unsigned int num_state_vars)
  : fEquation_Rhs(equation),
    fIntegrationOrder(integrationOrder),
    fNoIntegrationVariables(num_integration_vars),
    fNoStateVariables(num_state_vars)
{
}

GUVIntegrationStepper::~GUVIntegrationStepper()
{
}

void GUVIntegrationStepper::ComputeRightHandSide( const double y[], double charge, double dydx[] ) 
{
   this->RightHandSide( y, charge, dydx );
}

void GUVIntegrationStepper::SetEquationOfMotion(GUVEquationOfMotion* newEquation)
{
  if( newEquation != 0 )
  {
    fEquation_Rhs= newEquation;
  }
}
