//
//  Simple interface class to GUIntegrationDriver (with does Runge Kutta integration)
//   that follows the interface of TGeoHelix
//
#include <iostream>  // for  cout / cerr 

#include "GUFieldPropagator.h"
#include "GUVEquationOfMotion.h"
#include "GUVIntegrationStepper.h"
#include "GUIntegrationDriver.h"

// #include "TMagFieldEquation.h"
// #include "TClassicalRK4.h"

GUFieldPropagator::GUFieldPropagator(GUIntegrationDriver* driver) // (GUVField* field)
   : fInitialPosition(), fInitialDirection(),
     fMomentumMag(0.0), fInitialCurvature(0.0), fStepLength(0.0)
     //, fCurrentPoint(0.0, 0.0, 0.0), fCurrentDirection(0.0, 0.0, 0.0)
{
   // Must create the Driver, Stepper and Equation ??

   constexpr int NumEq= 6;
   // TMagFieldEquation *pEquation = TMagFieldEquation<ConstMagField,NumEq>(uniformField);
   // GUVEquationOfMotion*  pEquation= EquationFactory::CreateMagEquation(field, NumEq);
   // GUVIntegrationStepper = new TClassicalRK4<pEquation,NumEq>;
   // fDriver  = new GUIntegrationDriver();
   fDriver = driver;

   for( int i=0; i<3; i++) {
      fCurrentPoint[i]= 0.0;
      fCurrentDirection[i]= 0.0;
   }
}

// Make a step from current point along the path and compute new point, direction and angle
void GUFieldPropagator::Step(double step)
{
  // Do the work HERE
  fStepLength= step;

  GUFieldTrack fieldTr( fInitialPosition, 
                        fInitialDirection,
                        fMomentumMag,
                        // fRestMass,
                        fCharge, 
                        0.0,  // time
                        0.0); // s_0  

  // Call the driver HERE
  
  // fInitialCurvature; 
  ThreeVector endPoint= fieldTr.GetPosition();
  ThreeVector endDir  = fieldTr.GetMomentumDirection();

  for ( int i=0; i<3; i++) { 
     fCurrentPoint[i]=     endPoint[i];
     fCurrentDirection[i]= endDir[i];
  }

  std::cerr << " Incomplete METHOD -- must call Driver / etc for integratio. "  << std::endl;
  exit(1);
}

// static std::vector<GUFieldPropagator*> fFieldPropagatorVec;
// May change to c-array for CUDA ... but likely CPU only

/// --------------  GUFieldPropagatorPool ------------------------------------
// #include "GUFieldPropagatorPool.h"   // For now, not a separate file

// static
GUFieldPropagatorPool* 
GUFieldPropagatorPool::Instance()
{
   // A lock is REQUIRED for the next line - TODO
   static GUFieldPropagatorPool sInstance;

   return &sInstance;
}

#if 0
//// ---------------------  Postpone handling of multiple 
GUFieldPropagator* 
GUFieldPropagatorPool::CreateOrFind( int noNeeded ) // , void** banks )
{
  static int numberCreated= -1;
  static GUFieldPropagatorPool* pInstance= Instance();

  // A lock is REQUIRED for this section - TODO
  if( numberCreated < noNeeded)
  {
    Extend(noNeeded);
    assert( fFieldPropagatorVec.size() == noNeeded );
    // fNum = fFieldPropagatorVec.size();
    numberCreated= noNeeded;
  }
}

GUFieldPropagator* GUFieldPropagatorPool::GetPropagator(int num)
{
   assert(num>=0);
   assert(num< fFieldPropagatorVec.size());
  
   return fFieldPropagatorVec[num];
}
#endif

#if 0
void GUFieldPropagatorPool::Extend(int noNeeded)
{
    int num= fFieldPropagatorVec.size();
    while ( num < noNeeded )
    {
      //  if( (banks != 0) && (banks[num]!=0) )
      //  fFieldPropagatorVec.push( new(banks(num)) GUFieldPropagator() );
      //  else
      fFieldPropagatorVec.push_back( new GUFieldPropagator() );
    }
}
#endif 
