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
                                                     int num_state_vars)
  : fAbstrEquation(equation),
    fIntegrationOrder(integrationOrder),
    fNoIntegrationVariables(num_integration_vars),
    fNoStateVariables(num_state_vars > 0 ? num_state_vars : num_integration_vars)
{
}

GUVIntegrationStepper::~GUVIntegrationStepper()
{
}

// This allows the method to cache the value etc - Not needed for now
// void GUVIntegrationStepper::ComputeRightHandSide( const double y[], double charge, double dydx[] ) 
// {
//    this->RightHandSide( y, charge, dydx );
// }

void GUVIntegrationStepper::SetEquationOfMotion(GUVEquationOfMotion* newEquation)
{
  if( newEquation != 0 )
  {
    fAbstrEquation= newEquation;
  }
}
