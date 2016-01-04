//===----------------------------------------------------------------------===//
/**
 * @file   StepperFactory.h
 * @brief  Abstract field class for Geant-V prototype
 * @author John Apostolakis
 * @date   12 January 2016
 */
//===----------------------------------------------------------------------===//

#ifndef FIELDPROPAGATOR_FACTORY_H
#define FIELDPROPAGATOR_FACTORY_H 1

#include <ostream>

// #include "base/inc/Error.h"
#include "Geant/Error.h"

#include "GUFieldPropagator.h"
#include "GUFieldPropagatorPool.h"
#include "FieldEquationFactory.h"
#include "StepperFactory.h"
#include "GUIntegrationDriver.h"


// template<typename Field_t> // , typename Equation_t>
class FieldPropagatorFactory
{
public:
   
  static constexpr unsigned int  Nvar = 6; // Integration will occur over 3-position & 3-momentum coord.
  // using Equation_t = TMagFieldEquation<Field_t,Nvar>;

  // Initialise the classes required for tracking in field

  static constexpr double   fDefaultMinStep      = 0.0001;  // = 1 micron
  static constexpr double   fDefaultEpsTolerance = 1.0e-4;

  template<typename Field_t> // , typename Equation_t>
  static GUFieldPropagator* CreatePropagator(Field_t&    gvField,
                               double      relativeEpsTolerance,
                               double      minStepSize= fDefaultMinStep);
  // Create a propagator for RK integration of the motion in the Field 'gvField'
  // Then register it with the Pool (as the prototype)
  //
  // The Field_t object which is passed must be on the heap.
  //  It will be owned by the Propagator

  static GUFieldPropagator* CreatePropagator( // Field_t&    gvField,
                          //   Equation_t* gvEquation=  nullptr,
                               GUIntegrationDriver&   integrDriver,
                               double                 relTol= fDefaultEpsTolerance);
   // The GUIntegrationDriver object which is passed must be on the heap.
   //  It will be owned by the Propagator

private:
   static void RegisterPropagator(GUFieldPropagator*);
};

// template<typename Field_t> // , typename Equation_t>
GUFieldPropagator*
FieldPropagatorFactory::CreatePropagator( // Field_t&              gvField,
                                          GUIntegrationDriver&  integrDriver,
                                          double                relEpsilonTolerance)
{
  // using Equation_t =  TMagFieldEquation<Field_t,Nvar>;
  GUFieldPropagator* fieldPropagator = nullptr;

  // constexpr double epsTol = 3.0e-4;               // Relative error tolerance of integration
  // int  statisticsVerbosity= 0;

  // GUFieldPropagator *
  fieldPropagator =
     new GUFieldPropagator(&integrDriver, relEpsilonTolerance);  // epsTol);

  // cout << " - Integration constraint:  eps_tol= " << relEpsilonTolerance << endl;
  Geant::Print("FieldPropagatorFactory::CreatePropagator",  
               "Parameters for RK integration in magnetic field: \n - Integration constraint:  eps_tol=  %8.3g\n",
               relEpsilonTolerance); 
        
  RegisterPropagator(fieldPropagator);

  return fieldPropagator;
}

template<typename Field_t> // , typename Equation_t>
GUFieldPropagator*
FieldPropagatorFactory::CreatePropagator(Field_t& gvField,
                                         double   relEpsilonTolerance,
                                         double   minStepSize)
{
  using Equation_t =  TMagFieldEquation<Field_t,Nvar>;
  // const char* method="FieldPropagatorFactory::CreatePropagator";
  // cout << *method << " called. " << endl;
  auto gvEquation = 
     FieldEquationFactory::CreateMagEquation<Field_t>(&gvField);
  // cout << "Parameters for RK integration in magnetic field: "; //  << endl;
  // cout << " - Driver minimum step (h_min) = " << minStepSize << endl;

  Geant::Print("FieldPropagatorFactory::CreatePropagator",  
               // "Parameters for RK integration in magnetic field: "
               " - Driver minimum step (h_min) = %8.3g\n",
               minStepSize); 
  
  auto // GUVIntegrationStepper*
     aStepper = StepperFactory::CreateStepper<Equation_t>(gvEquation); // Default stepper

  int   statisticsVerbosity= 0;
  auto integrDriver = new GUIntegrationDriver( minStepSize,
                                               aStepper,
                                               Nvar,
                                               statisticsVerbosity);

  return CreatePropagator( *integrDriver, relEpsilonTolerance );
}

// template<typename Field_t, typename Equation_t>
void
FieldPropagatorFactory::RegisterPropagator(GUFieldPropagator* fieldPropagator)
{
  GUFieldPropagatorPool* fpPool= GUFieldPropagatorPool::Instance();
  assert( fpPool );  // Cannot be zero
  if( fpPool ) {
     fpPool->RegisterPrototype( fieldPropagator );
     // printf( "FieldPropagatorFactory: Registered Prototype field-prop %p\n", fieldPropagator );
  } else {
     Geant::Error("PrepareRkIntegration","Cannot find GUFieldPropagatorPool Instance.");
  }
}
#endif
