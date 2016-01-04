//
//  Simple interface class to GUIntegrationDriver (with does Runge Kutta integration)
//   that follows the interface of TGeoHelix
//
#include <iostream>  // for  cout / cerr 

#include "GUFieldPropagator.h"

// #include "GUVEquationOfMotion.h"
#include "TMagFieldEquation.h"
#include "GUVIntegrationStepper.h"
#include "GUIntegrationDriver.h"
#include "GUVEquationOfMotion.h"

#include "TMagFieldEquation.h"
#include "TClassicalRK4.h"

using ThreeVector = vecgeom::Vector3D<double>;

GUFieldPropagator::GUFieldPropagator(GUIntegrationDriver* driver, double eps)
  : fDriver(driver), fEpsilon(eps)
{
}

// ToDo-s/ideas:
//  - Factory to create the Driver, Stepper and Equation

template<typename FieldType>  // , typename StepperType>
GUFieldPropagator::GUFieldPropagator(FieldType* magField, double eps, double hminimum)
   : fEpsilon(eps)
{
   constexpr int NumEq= 6;
   using  EquationType=  TMagFieldEquation<FieldType, NumEq>;
   
   int statVerbose= 1;
   auto *pEquation = new EquationType(magField, NumEq);
      // new TMagFieldEquation<FieldType,NumEq>(magField, NumEq);

   // auto stepper = new StepperType<GvEquationType,NumEq>(gvEquation);
   auto stepper =      new TClassicalRK4<EquationType,NumEq>(pEquation);      
   auto integrDriver = new GUIntegrationDriver( hminimum,
                                               stepper,
                                               NumEq,
                                               statVerbose);
   fDriver= integrDriver;
}

GUFieldPropagator* GUFieldPropagator::Clone() const 
{
   return new GUFieldPropagator( fDriver->Clone(), fEpsilon );
}

// Make a step from current point along the path and compute new point, direction and angle
// VECCORE_ATT_HOST_DEVICE                 
bool
GUFieldPropagator::DoStep( ThreeVector const & startPosition, ThreeVector const & startDirection,
                                   int const & charge,             double const & startMomentumMag,
                                double const & step,
                           ThreeVector       & endPosition,
                           ThreeVector       & endDirection
         )
{
  // Do the work HERE
  GUFieldTrack yTrackIn( startPosition, 
                        startDirection * startMomentumMag,
                        // fCharge, 
                        0.0); // s_0  xo
  GUFieldTrack yTrackOut( yTrackIn );
  
  // Call the driver HERE
  fDriver->InitializeCharge( charge );
  bool goodAdvance=
     fDriver->AccurateAdvance( yTrackIn, step, fEpsilon, yTrackOut ); // , hInitial );

  // fInitialCurvature; 
  endPosition=  yTrackOut.GetPosition();
  endDirection= yTrackOut.GetMomentumDirection();
  return goodAdvance;
}

GUVField* GUFieldPropagator::GetField() 
{
   GUVField* pField = nullptr;
   auto driver= GetIntegrationDriver();
   if( driver ){
     auto equation= driver->GetEquationOfMotion();
     if( equation ) {
       pField = equation->GetFieldObj();
     }
   }
   return pField;
}

// static std::vector<GUFieldPropagator*> fFieldPropagatorVec;
// May change to c-array for CUDA ... but likely CPU only

