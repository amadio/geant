//===----------------------------------------------------------------------===//
/**
 * @file   FieldPropagatorFactory.h
 * @brief  Class to create Field Propagator objects for Geant-V prototype
 * @author John Apostolakis
 * @date   12 January 2016
 */
//===----------------------------------------------------------------------===//

#ifndef FIELDPROPAGATOR_FACTORY_H
#define FIELDPROPAGATOR_FACTORY_H 1

#include <ostream>

// #include "base/inc/Geant/Error.h"
// #include "Geant/Error.h"

#include "GUFieldPropagator.h"
#include "GUFieldPropagatorPool.h"

// #ifndef FLEXIBLE_FIELD
#include "ScalarMagFieldEquation.h"
#include "FieldEquationFactory.h"
#include "StepperFactory.h"
#include "ScalarIntegrationDriver.h"
// #else
#include "MagFieldEquation.h"
#include "CashKarp.h"
#include "FlexIntegrationDriver.h"
#include "SimpleIntegrationDriver.h"
// #endif

// template<typename Field_t> // , typename Equation_t>
class FieldPropagatorFactory
{
public:
   
  static constexpr unsigned int  Nvar = 6; // Integration will occur over 3-position & 3-momentum coord.
  // using Equation_t = TMagFieldEquation<Field_t,Nvar>;

  // Initialise the classes required for tracking in field

  static constexpr double   fDefaultMinStep      = 0.0001;  // = 1 micron
  static constexpr double   fDefaultEpsTolerance = 1.0e-4;

  /** brief Create a 'full' propagator with both scalar & vector/flexible drivers */
  template<typename Field_t>
  static
     GUFieldPropagator* CreatePropagator(Field_t&    gvField,
                                         double      relativeTolerance,
                                         double      minStep= fDefaultMinStep);
  // To be used for RK integration of the motion in the Field 'gvField'
  // Will register it with the Pool (as the prototype)
  // The Field_t object which is passed must be on the heap.
  //  ( It will be owned by the Propagator. )

  /** @ brief  Create object using given drivers (on the heap) - obtains their ownership. */
  static
    GUFieldPropagator* CreatePropagator( ScalarIntegrationDriver& integrDriver,
                                         double                   relTolerance,
                                         FlexIntegrationDriver*   flexDrv = nullptr
        );
   // The ScalarIntegrationDriver object which is passed must be on the heap.
   //  It will be owned by the Propagator

  /** @brief Create a 'scalar' propagator (with only scalar driver) for RK integration  */
  template<typename Field_t>  
  static GUFieldPropagator* CreateScalarPropagator(Field_t&    gvField,
                                                   double      relativeEpsTolerance,
                                                   double      minStepSize= fDefaultMinStep);
  // Registers it with the Pool (as the prototype)

private:
  //  Helper methods
  
  /** @ brief Auxiliary methods to create scalar driver */
  template<typename Field_t>
  static
    ScalarIntegrationDriver* CreateScalarDriver(Field_t&    gvField,
                                                double      relativeEpsTolerance,
                                                double      minStepSize= fDefaultMinStep);
  // Create a 'scalar' driver for RK integration of the motion in the Field 'gvField'.


  /** @ brief Auxiliary methods to create scalar driver */
  template<typename Field_t>
  static
    FlexIntegrationDriver* CreateFlexibleDriver(Field_t&    gvField,
                                                double      relativeEpsTolerance,
                                                double      minStepSize= fDefaultMinStep);
  
  
private:
   static void RegisterPropagator(GUFieldPropagator*);
};

//______________________________________________________________________________
// template<typename Field_t> // , typename Equation_t>
inline
GUFieldPropagator*
FieldPropagatorFactory::CreatePropagator( // Field_t&              gvField,
                                          ScalarIntegrationDriver&  integrDriver,
                                          double                    relEpsilonTolerance,
                                          FlexIntegrationDriver*    flexDriver
   )
{
  // using Equation_t =  TMagFieldEquation<Field_t,Nvar>;
  GUFieldPropagator* fieldPropagator = nullptr;

  // constexpr double epsTol = 3.0e-4;               // Relative error tolerance of integration

  // GUFieldPropagator *
  fieldPropagator =
     new GUFieldPropagator(&integrDriver, relEpsilonTolerance, flexDriver );

  std::cout << " - Integration constraint:  eps_tol= " << relEpsilonTolerance << std::endl;
  // Geant::Printf("FieldPropagatorFactory::CreatePropagator",  
  //             "Parameters for RK integration in magnetic field: \n - Integration constraint:  eps_tol=  %8.3g\n",
  //              relEpsilonTolerance); 
        
  RegisterPropagator(fieldPropagator);

  return fieldPropagator;
}

//______________________________________________________________________________
template<typename Field_t>
inline
GUFieldPropagator*
FieldPropagatorFactory::CreatePropagator(Field_t&  gvField,
                                         double    relativeTolerance,
                                         double    minStep )
{
   auto  // ScalarIntegrationDriver
      scalarDriver = CreateScalarDriver( gvField, relativeTolerance, minStep );
   FlexIntegrationDriver*  // auto
      flexibleDriver = CreateFlexibleDriver( gvField, relativeTolerance, minStep );
   return FieldPropagatorFactory::CreatePropagator( *scalarDriver, relativeTolerance, flexibleDriver );
}


//______________________________________________________________________________
template<typename Field_t> // , typename Equation_t>
inline
ScalarIntegrationDriver*
FieldPropagatorFactory::CreateScalarDriver(Field_t& gvField,
                                           double   /*relEpsilonTolerance*/,
                                           double   minStepSize)
{
  const char* methodName="FieldPropagatorFactory::CreateScalarDriver";
  int   statisticsVerbosity= 0;
  
  // cout << methodName << " called. " << endl;

  using Equation_t =  ScalarMagFieldEquation<Field_t, Nvar>;
  auto gvEquation = 
     FieldEquationFactory::CreateMagEquation<Field_t>(&gvField);
  auto // VScalarIntegrationStepper*
     aStepper = StepperFactory::CreateStepper<Equation_t>(gvEquation); // Default stepper

  auto scalarDriver = new ScalarIntegrationDriver( minStepSize,
                                                   aStepper,
                                                   Nvar,
                                                   statisticsVerbosity);
  
  std::cout << methodName << ": Parameters for RK integration in magnetic field: "; //  << endl;
  std::cout << " - Driver minimum step (h_min) = " << minStepSize << std::endl;

  // Geant::Print(methodName, 
               // "Parameters for RK integration in magnetic field: "
  //            " - Driver minimum step (h_min) = %8.3g\n", minStepSize); 
  
  return scalarDriver;
}

//______________________________________________________________________________
template<typename Field_t>
inline
FlexIntegrationDriver*
FieldPropagatorFactory::CreateFlexibleDriver(Field_t&    gvField,
                                            double   /*relEpsilonTolerance*/,
                                            double   minStepSize)
{
  const char* methodName="FieldPropagatorFactory::CreateFlexibleDriver";
  int   statsVerbose= 1;
  
  // cout << methodName << " called. " << endl;

  // New flexible (scalar + vector) versions of field, equation, ...
  constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec
  
  using Equation_t =  MagFieldEquation<Field_t>;  // Flexible version  
  auto gvEquation = new Equation_t(&gvField);

  using StepperType = CashKarp<Equation_t, Nposmom>;
  auto myStepper = new StepperType(gvEquation);
  // new CashKarp<GvEquationType,Nposmom>(gvEquation);

  using  DriverType = SimpleIntegrationDriver<StepperType,Nposmom>;
  auto vectorDriver = new DriverType( minStepSize,
                                      myStepper,
                                      Nposmom,
                                      statsVerbose);

  std::cout << methodName << ": Parameters for RK integration in magnetic field: "
            << " - Driver minimum step (h_min) = " << minStepSize << std::endl;

  // Geant::Print(methodName,
               // "Parameters for RK integration in magnetic field: "
  //            " - Driver minimum step (h_min) = %8.3g\n", minStepSize);
               
  
  return vectorDriver;
}



// template<typename Field_t, typename Equation_t>

#endif
