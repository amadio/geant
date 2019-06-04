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

#include "Geant/GUFieldPropagator.h"
#include "Geant/GUFieldPropagatorPool.h"

// #ifndef FLEXIBLE_FIELD
#include "Geant/ScalarMagFieldEquation.h"
#include "Geant/FieldEquationFactory.h"
#include "Geant/StepperFactory.h"
#include "Geant/ScalarIntegrationDriver.h"
// #else
#include "Geant/MagFieldEquation.h"
#include "Geant/CashKarp.h"
#include "Geant/DormandPrince5RK.h"
#include "Geant/BogackiShampine23RK.h"

#include "Geant/FlexIntegrationDriver.h"
#include "Geant/SimpleIntegrationDriver.h"

#include "Geant/OldIntegrationDriver.h"
#include "Geant/RollingIntegrationDriver.h"
// #endif

enum StepperTypeNum { kBogackiShampineStepper=4 ,
                      kCashKarpStepper= 6,
                      kDormandPrince45Stepper= 7,
                      // Dormand Prince with fixed type of Driver:
                      kDormandPrince45StepperOldDriver= 8,
                      kDormandPrince45StepperSimpleDriver= 9,
                      kDormandPrince45StepperRollingDriver= 10,
                      // Auxiliary values
                      kMinStepperOrder = 3,                      
                      kMaxStepperOrder = 11,
                      kUndefinedStepperType = 0 };

// template<typename Field_t> // , typename Equation_t>
class FieldPropagatorFactory {
public:
  static constexpr unsigned int NumVar = 6; // Integration will occur over 3-position & 3-momentum coord.
  // using Equation_t = TMagFieldEquation<Field_t,NumVar>;

  // Initialise the classes required for tracking in field

  static constexpr double fDefaultMinStep      = 0.0001; // = 1 micron
  static constexpr double fDefaultEpsTolerance = 1.0e-4;

  /** brief Create a 'full' propagator with both scalar & vector/flexible drivers */
  template <typename Field_t>
  static GUFieldPropagator *CreatePropagator(Field_t & gvField,
                                             double    relativeTolerance,
                                             double    minStep, // = fDefaultMinStep,
                                             int       stepperTypeNo /*= 0*/ );

  // To be used for RK integration of the motion in the Field 'gvField'
  // Will register it with the Pool (as the prototype)
  // The Field_t object which is passed must be on the heap.
  //  ( It will be owned by the Propagator. )

  /** @ brief  Create object using given drivers (on the heap) - obtains their ownership. */
  static GUFieldPropagator *CreatePropagator(ScalarIntegrationDriver *integrDriver, double relTolerance,
                                             FlexIntegrationDriver *flexDriver = nullptr);
  // The ScalarIntegrationDriver object which is passed must be on the heap.
  //  It will be owned by the Propagator

  /** @brief Create a 'scalar' propagator (with only scalar driver) for RK integration  */
  template <typename Field_t>
  static GUFieldPropagator *CreateScalarPropagator(Field_t &gvField,
                                                   double relativeEpsTolerance,
                                                   double minStepSize ); // = fDefaultMinStep);
  // Registers it with the Pool (as the prototype)

private:
  //  Helper methods
  /** @  brief Obtain enum of stepper type - or default */
  static StepperTypeNum GetStepperTypeId( int stepperTypeNo );
  
  /** @ brief Auxiliary methods to create scalar driver */
  template <typename Field_t>
  static ScalarIntegrationDriver *CreateScalarDriver(Field_t & gvField,
                                                     double    relativeEpsTolerance,
                                                     double    minStepSize ); // = fDefaultMinStep);
  // Create a 'scalar' driver for RK integration of the motion in the Field 'gvField'.

  /** @ brief Auxiliary methods to create flexible driver */
  template <typename Field_t>
  static FlexIntegrationDriver *CreateFlexibleDriver(Field_t &gvField,
                                                     double relativeEpsTolerance,
                                                     double minStepSize, // = fDefaultMinStep,
                                                     StepperTypeNum stepperTypeId // = kCashKarpStepper
                                                    );

  /** @ brief Auxiliary methods to create flexible driver */  
template <typename Equation_t, typename Stepper_t>
  static FlexIntegrationDriver*  CreateDriverForStepper(Equation_t & equation,
                                                        double minStepSize,
                                                        double relEpsilonTolerance,                                                                 
                                                        Stepper_t **stepperObj // = nullptr
                                                       );

public:
  static bool fVerboseConstruct;
  // Verbosity for construction

private:
  static void RegisterPropagator(GUFieldPropagator *);
};

//______________________________________________________________________________
// template<typename Field_t> // , typename Equation_t>
inline GUFieldPropagator *FieldPropagatorFactory::
    CreatePropagator( // Field_t&               gvField,
                      ScalarIntegrationDriver * integrDriver,
                      double                    relEpsilonTolerance,
                      FlexIntegrationDriver   * flexDriver )
{
  // using Equation_t =  TMagFieldEquation<Field_t,NumVar>;
  const char *methodName             = "FieldPropagatorFactory::CreatePropagator";
  GUFieldPropagator *fieldPropagator = nullptr;

  // constexpr double epsTol = 3.0e-4;               // Relative error tolerance of integration

  assert(integrDriver); // Cannot be null!
  if (fVerboseConstruct)
    std::cout << "Check scalar Driver: max Num steps= " << integrDriver->GetMaxNoSteps() << std::endl;

  // GUFieldPropagator *
  fieldPropagator = new GUFieldPropagator(integrDriver, relEpsilonTolerance, flexDriver);

  assert( flexDriver->GetMaxRelativeEpsilon() == relEpsilonTolerance );
  if( flexDriver->GetMaxRelativeEpsilon() != relEpsilonTolerance ) {
     std::cerr << "ERROR> Difference in Epsilon value requested ( = " << relEpsilonTolerance << " ) " 
        " , and value stored in Driver ( value = " << flexDriver->GetMaxRelativeEpsilon() << std::endl;
     exit(1);     
  }
  
  if (fVerboseConstruct) {
    std::cout << methodName << " ( scalar, double, flex-driver ) called "
              << " - Integration constraint:  eps_tol= " << relEpsilonTolerance << std::endl;
    std::cout << methodName << "  scalarDriver = " << &integrDriver << std::endl;
    std::cout << methodName << "  vectorDriver = " << flexDriver << std::endl;
    // geant::Printf("FieldPropagatorFactory::CreatePropagator",
    //             "Parameters for RK integration in magnetic field: \n - Integration constraint:  eps_tol=  %8.3g\n",
    //              relEpsilonTolerance);
  }

  RegisterPropagator(fieldPropagator);

  return fieldPropagator;
}

//______________________________________________________________________________
template <typename Field_t>
inline GUFieldPropagator *FieldPropagatorFactory::CreatePropagator(Field_t &gvField,
                                                                   double relativeTolerance,
                                                                   double minStep,
                                                                   int    stepperTypeNo
   )
{
  const char *methodName = "FieldPropagatorFactory::CreatePropagator";
  const char *methodSig  = "( templated<Field_t> field, double, double )";
  auto stepperTypeId= GetStepperTypeId( stepperTypeNo );

  if (fVerboseConstruct)
     std::cout << methodName << " " << methodSig << " called with stepperId="
               << stepperTypeNo << std::endl;

  auto scalarDriver = CreateScalarDriver(gvField, relativeTolerance, minStep); 
      
  // FlexIntegrationDriver * 
  auto flexibleDriver = CreateFlexibleDriver(gvField, relativeTolerance, minStep, stepperTypeId);
                              
  return FieldPropagatorFactory::CreatePropagator(scalarDriver, relativeTolerance, flexibleDriver);
}

StepperTypeNum
inline
FieldPropagatorFactory::GetStepperTypeId( int stepperTypeNo )
{
  StepperTypeNum stepperTypeId= kUndefinedStepperType;
  if( stepperTypeNo == 3 || stepperTypeNo == 4 ) stepperTypeId = kBogackiShampineStepper;
  else if ( stepperTypeNo == 5 )  stepperTypeId =  kCashKarpStepper;  
  else if ( stepperTypeNo == 7 )  stepperTypeId = kDormandPrince45Stepper;
  else {
     std::cerr << "FieldPropagationFactory: Default type of Stepper chosen - as none of the expected values requested." << std::endl;
     stepperTypeId =  kCashKarpStepper;
  }
  return stepperTypeId;
}

//______________________________________________________________________________
template <typename Field_t> // , typename Equation_t>
inline ScalarIntegrationDriver *FieldPropagatorFactory::CreateScalarDriver(Field_t &gvField,
                                                                           double relEpsilonTolerance,
                                                                           double minStepSize)
{
  const char *methodName  = "FieldPropagatorFactory::CreateScalarDriver";
  int statisticsVerbosity = 0;

  // cout << methodName << " called. " << endl;

  using Equation_t = ScalarMagFieldEquation<Field_t, NumVar>;
  auto gvEquation  = FieldEquationFactory::CreateMagEquation<Field_t>(&gvField);
  auto                                                                  // VScalarIntegrationStepper*
      aStepper = StepperFactory::CreateStepper<Equation_t>(gvEquation); // Default stepper

  auto scalarDriver = new ScalarIntegrationDriver(minStepSize, aStepper, relEpsilonTolerance, NumVar, statisticsVerbosity);

  if (fVerboseConstruct) {
    std::cout << methodName << ": Parameters for RK integration in magnetic field: "; //  << endl;
    std::cout << " - Driver minimum step (h_min) = " << minStepSize << scalarDriver->GetMaxNoSteps() << std::endl;
    // Test the object ...

    // geant::Print(methodName,
    // "Parameters for RK integration in magnetic field: "
    //            " - Driver minimum step (h_min) = %8.3g\n", minStepSize);
  }

  return scalarDriver;
}


//______________________________________________________________________________
template <typename Equation_t, typename Stepper_t>
  inline FlexIntegrationDriver* FieldPropagatorFactory::
   CreateDriverForStepper(Equation_t&  equation,
                          double       relEpsilonTolerance,                          
                          double       minStepSize,
                          Stepper_t  **stepperObj )
{
  const char *methodName = "FieldPropagatorFactory::CreateDriverForStepper";

  // std::cout << methodName << " called. " << std::endl;

  // New flexible (scalar + vector) versions of field, equation, ...
  constexpr unsigned int Nposmom = 6; // Position 3-vec + Momentum 3-vec
  auto myStepper    = new Stepper_t(&equation);

  using DriverType  = SimpleIntegrationDriver<Stepper_t, Nposmom>;
  // using DriverType = OldIntegrationDriver<StepperType, Nposmom>;
  // using DriverType = RollingIntegrationDriver<StepperType, Nposmom>;
  
  auto vectorDriver = new DriverType(minStepSize, myStepper, relEpsilonTolerance, Nposmom);

  assert(vectorDriver);

  if (fVerboseConstruct) {
    std::cout << methodName << ": Parameters for RK integration in magnetic field: "
              << " - Driver minimum step (h_min) = " << minStepSize << std::endl;
    std::cout << methodName << ": created vector driver = " << vectorDriver << std::endl;
    // geant::Print(methodName,
    //              "Parameters for RK integration in magnetic field: "
    //             " - Driver minimum step (h_min) = %8.3g\n", minStepSize);
  }
  if( stepperObj ) *stepperObj= myStepper;
  
  return vectorDriver;
}

//______________________________________________________________________________
template <typename Field_t>
inline FlexIntegrationDriver* FieldPropagatorFactory::
   CreateFlexibleDriver(Field_t &gvField,
                        double relEpsilonTolerance,
                        double minStepSize,
                        StepperTypeNum stepperTypeId
   )
{
  using Equation_t = MagFieldEquation<Field_t>; // Flexible version
  constexpr unsigned int Nposmom = 6; // Position 3-vec + Momentum 3-vec
  using StepperTypeCK456    = CashKarp<Equation_t, Nposmom>;
  using StepperTypeDoPri457 = DormandPrince5RK<Equation_t, Nposmom>;
  using std::cout;
  using std::endl;
  
  FlexIntegrationDriver* vectorDriver= nullptr;
  
  const char *methodName = "FieldPropagatorFactory::CreateFlexibleDriver";

  std::cout << methodName << " - called with stepperTypeID = " << stepperTypeId;
  stepperTypeId = kDormandPrince45StepperRollingDriver;

  std::cout << " -- overrode to use " << stepperTypeId 
            << " ( = " << kDormandPrince45StepperRollingDriver << " ) "
            << " ,  i.e.  Rolling Driver with Dormand Prince 4/5 " 
            << std::endl;
  assert( stepperTypeId == kDormandPrince45StepperRollingDriver );
  
  // New flexible (scalar + vector) versions of field, equation, ...

  auto gvEquation  = new Equation_t(&gvField);

  StepperTypeDoPri457* myStepperDoPri5 = nullptr;
  StepperTypeCK456*    myStepperCK5    = nullptr;

  // using StepperTypeBS234 = BogackiShampine23RK<Equation_t, Nposmom>;  
  // StepperTypeBS234*    myStepperBS234  = nullptr;

  switch ( stepperTypeId )
  {
     case kDormandPrince45Stepper:
     case kMaxStepperOrder:
        vectorDriver = CreateDriverForStepper<Equation_t, StepperTypeDoPri457>
           (*gvEquation,  minStepSize, relEpsilonTolerance, &myStepperDoPri5);
        cout << methodName << ": created Dormand Prince 4/5 (7 stage) stepper with default driver." << endl;
        break;
/****
     case kMinStepperOrder:
     case kBogackiShampineStepper:
        myStepperBS234 = new StepperTypeBS234(gvEquation);
        vectorDriver = new SimpleIntegrationDriver<StepperTypeBS234, Nposmom>
           (minStepSize, myStepperBS234, relEpsilonTolerance, Nposmom );
        cout << methodName << ": created Dormand Prince 4/5 (7 stage) stepper with Simple driver. " << endl;
        break;
*****/
     case kCashKarpStepper:      
        myStepperCK5 = new StepperTypeCK456(gvEquation);
        vectorDriver = new SimpleIntegrationDriver<StepperTypeCK456, Nposmom>
                                     (minStepSize, myStepperCK5, relEpsilonTolerance, Nposmom );
        cout << methodName << ": created Cash Karp 4/5 (7 stage) stepper with Simple driver." << endl;        
        break;
/******        
       // Explicit request for alternative *Driver* type:  Old Driver (earlier iteration of Simple)
     case kDormandPrince45StepperOldDriver:
           myStepperDoPri5 = new StepperTypeDoPri457(gvEquation);
           vectorDriver = new OldIntegrationDriver<StepperTypeDoPri457, Nposmom>
                                  (minStepSize, myStepperDoPri5, relEpsilonTolerance, Nposmom);
           cout << methodName << ": created Dormand Prince 4/5 (7 stage) stepper with Old Simple driver." << endl;
        break;
****/
       // Explicit request for alternative *Driver* type:   'Simple' Driver
     case kDormandPrince45StepperSimpleDriver:
           myStepperDoPri5 = new StepperTypeDoPri457(gvEquation);
           vectorDriver = new SimpleIntegrationDriver<StepperTypeDoPri457, Nposmom>
                                  (minStepSize, myStepperDoPri5, relEpsilonTolerance, Nposmom);
           cout << methodName << ": created Dormand Prince 4/5 (7 stage) stepper with     Simple driver." << endl;           
        break;

       // Explicit request for alternative *Driver* type:   'Simple' Driver
     case kDormandPrince45StepperRollingDriver:
           myStepperDoPri5 = new StepperTypeDoPri457(gvEquation);
           vectorDriver = new RollingIntegrationDriver<StepperTypeDoPri457, Nposmom>
                                  (minStepSize, myStepperDoPri5, relEpsilonTolerance, Nposmom);
           cout << methodName << ": created Dormand Prince 4/5 (7 stage) stepper with    Rolling driver." << endl;                      
        break;

    case kUndefinedStepperType:
        std::cerr << methodName << " FATAL ERROR : Stepper type is explitictly undefined." << std::endl;
       // Old: "- using DormandPrince (default.)" << std::endl;
        exit(1);
        
     default:
        std::cerr << methodName << " WARNING : Stepper type is another value - case 'default' - using CashKarp." << std::endl;
        myStepperCK5 = new StepperTypeCK456(gvEquation);
        vectorDriver = new SimpleIntegrationDriver<StepperTypeCK456, Nposmom>
                                     (minStepSize, myStepperCK5, relEpsilonTolerance, Nposmom );
        break;        
        
  }
  // using DriverType  = SimpleIntegrationDriver<StepperType, Nposmom>;
  // auto vectorDriver = new DriverType(minStepSize, myStepper, Nposmom);

  assert(vectorDriver);

  if (fVerboseConstruct) {
    std::cout << methodName << ": Parameters for RK integration in magnetic field: "
              << " - Driver minimum step (h_min) = " << minStepSize << std::endl;
    std::cout << methodName << ": created vector driver = " << vectorDriver << std::endl;
    // geant::Print(methodName,
    //              "Parameters for RK integration in magnetic field: "
    //             " - Driver minimum step (h_min) = %8.3g\n", minStepSize);
  }

  return vectorDriver;
}

// template<typename Field_t, typename Equation_t>

#endif
