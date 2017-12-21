//
//  Simple interface class to ScalarIntegrationDriver (with does Runge Kutta integration)
//   that follows the interface of TGeoHelix
//
#include <iostream>  // for  cout / cerr 

#include "GUFieldPropagator.h"

#include "VScalarEquationOfMotion.h"
#include "ScalarMagFieldEquation.h"
#include "VScalarIntegrationStepper.h"
#include "ScalarIntegrationDriver.h"

#include "GUTCashKarpRKF45.h"  //  ie ScalarCashKarp

#include "MagFieldEquation.h"
#include "CashKarp.h"
#include "FlexIntegrationDriver.h"
#include "SimpleIntegrationDriver.h"

// template <class Equation, unsigned int> using ScalarCashKarp= GUTCashKarpRKF45;

using ThreeVector = vecgeom::Vector3D<double>;

FlexIntegrationDriver*  GUFieldPropagator::fVectorDriver= nullptr;
double                  GUFieldPropagator::fEpsilon = 1.0e-4;
//------------------------------------------------------------------------------------
GUFieldPropagator::GUFieldPropagator(ScalarIntegrationDriver* driver,
                                     double eps,
                                     FlexIntegrationDriver* flexDriver )
  : fScalarDriver(driver),
    fVerboseConstruct(false)
    // fEpsilon(eps)
{
   const char* methodName="GUFieldPropagator constructor (scalarDriver, eps, *flex - with default = null)";

   fEpsilon= eps;


   if( fVerboseConstruct ) {
      std::cout << "------------------------------------------------------------------------------------" << std::endl;
      std::cout << methodName << std::endl
                << "  Arguments:  scalarDrv= " << driver << " eps = " << eps << " flexDriver = " << flexDriver
                << std::endl;
      std::cout << "  Existing vectorDriver (static) = " << fVectorDriver << std::endl;
   }
   
   if ( !fVectorDriver && flexDriver ) {      
      SetFlexIntegrationDriver( flexDriver );
      if( fVerboseConstruct ) {
         std::cout << " Installing flex/vectorDriver = " << flexDriver << std::endl;
         std::cout << " Confirmed  address           = " << GetFlexibleIntegrationDriver() << std::endl;
      }
      assert( flexDriver == GetFlexibleIntegrationDriver() ); 
   } else {
      if ( fVectorDriver && flexDriver )
         std::cout << "GUFieldPropagator> Not overwriting Vector/Flexible Driver" << std::endl;
      if ( !fVectorDriver && !flexDriver ) {               
         std::cout << "==============================================================================" << std::endl;
         std::cout << "WARNING: Gufieldpropagator> Not setting Vector/Flexible Driver: none provided." << std::endl;
         std::cout << "==============================================================================" << std::endl;
      }
   }

   if( fVerboseConstruct ) {
      std::cout << "GUFieldPropagator constructor> ptr= " << this << std::endl
                << " scalar driver = " << fScalarDriver
                << " flex/vector driver = " << flexDriver << std::endl;
   }
}

//------------------------------------------------------------------------------------
GUFieldPropagator::~GUFieldPropagator()
{
   if( fVerboseConstruct )
      std::cout << "GUFieldPropagator destructor called for ptr= " << this << std::endl;
}

//------------------------------------------------------------------------------------
void GUFieldPropagator::SetFlexIntegrationDriver( FlexIntegrationDriver * flexDriver)
{
   const std::string methodName = "GUFieldPropagator::SetFlexIntegrationDriver";
   if ( !fVectorDriver && flexDriver ) {
      std::cout << "Replacing Vector/Flexible Driver" << std::endl;
      fVectorDriver = flexDriver;
   }
   else if ( !flexDriver )
   {
      std::cout << methodName // << "GUFieldPropagator::SetFlexIntegrationDriver",
         << "> Not overwriting with Vector/Flexible Driver" << std::endl;
   }
}

// ToDo-s/ideas:
//  - Factory to create the Driver, Stepper and Equation

//____________________________________________________________________________________
template<typename FieldType>  // , typename StepperType>
GUFieldPropagator::GUFieldPropagator(FieldType* magField, double eps, double hminimum)
  // : fEpsilon(eps)
{
   constexpr unsigned int Nposmom = 6; // Number of Integration variables - 3 position, 3 momentum

   if( 0.0 < eps && eps < 1.0 ) 
      fEpsilon= eps;
   
// #if 0
   using  ScalarEquationType=  ScalarMagFieldEquation<FieldType, Nposmom>;
   int statVerbose= 1;
   auto *pEquation = new ScalarEquationType(magField, Nposmom);
      // new ScalarFieldEquation<FieldType,Nposmom>(magField, Nposmom);

   // auto stepper = new StepperType<ScalarEquationType,Nposmom>(gvEquation);
   auto scalarStepper =      new // ScalarCashKarp
           GUTCashKarpRKF45<ScalarEquationType,Nposmom>(pEquation);
   auto scalarDriver = new ScalarIntegrationDriver( hminimum,
                                                    scalarStepper,
                                                    Nposmom,
                                                    statVerbose);
   fScalarDriver= scalarDriver;
// #else
//    fScalarDriver= nullptr;
// #endif

   if( ! fVectorDriver ) {    // Share it between all instances
      
      // Create the flexible (vector or scalar) objects 
      using FlexEquationType = MagFieldEquation<FieldType>;
      auto  gvEquation = new FlexEquationType(magField);
      using FlexStepperType = CashKarp<FlexEquationType,Nposmom>;
      auto myFlexStepper = new FlexStepperType(gvEquation);
      int statsVerbose=1;
      auto flexDriver =
         new SimpleIntegrationDriver<FlexStepperType,Nposmom> (hminimum,
                                                               myFlexStepper,
                                                               Nposmom,
                                                               statsVerbose);
      fVectorDriver= flexDriver;
   }

     std::cout << "GUFieldPropagator constructor > this = " << this << std::endl;
}

// #ifdef FP_CLONE_METHOD
GUFieldPropagator* GUFieldPropagator::Clone() const 
{
    return new GUFieldPropagator( fScalarDriver->Clone(),
                                  fEpsilon );
}
// #endif

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
  const char *methodName ="GUFieldPropagator::DoStep";
  bool goodAdvance=false;
  bool verbose= false;

  verbose= true;
  static bool infoPrinted= false;
  
// #define USE_FLEXIBLE_FOR_SCALAR 1
  
#ifdef USE_FLEXIBLE_FOR_SCALAR
  // Try using the flexible driver first - for a single value ... ( to be improved )
  double     vals[6] = { startPosition.x(), startPosition.y(), startPosition.z(),
                         startMomentumMag * startDirection.x(),
                         startMomentumMag * startDirection.y(),
                         startMomentumMag * startDirection.z() };
  FieldTrack yTrackInFT( vals ), yTrackOutFT;
  double invMomentumMag = 1.0 / startMomentumMag;
  double chargeFlt= charge;
  bool   okFlex= false;

  // Checks ...
  assert( fVectorDriver && "Vector Driver is not SET" );
  if( verbose ) { 
     std::cout << methodName << " Dump of arguments: " << std::endl;
     std::cout << "   Step       = " << step << std::endl;
     std::cout << "   yTrackInFT = " << yTrackInFT << std::endl;
  }


#ifdef EXTEND_SINGLE
  if( verbose && !infoPrinted ) {
     std::cout << methodName << " > Using Flexible/Vector Driver - for 1 track " << fScalarDriver << std::endl;
     infoPrinted = true;
  }
  
  // Using the capabiity of the flexible driver to integrate a single track -- New 25.01.2018
  fVectorDriver->AccurateAdvance( yTrackInFT, step, chargeFlt, fEpsilon, yTrackOutFT, okFlex );
#else
  if( verbose && !infoPrinted ) {
     std::cout << methodName << " > Using VectorDriver ( Real_v ) with 1 track" << fScalarDriver << std::endl;
     infoPrinted = true;
  }

  // Harnessing the vector driver to integrate just a single track -- MIS-USE
  fVectorDriver->AccurateAdvance( &yTrackInFT, &step, &chargeFlt, fEpsilon, &yTrackOutFT, 1, &okFlex );
#endif
  if( verbose ) 
     std::cout << " Results:  good = " << okFlex << " track out= " << yTrackOutFT << std::endl;
// #endif

  // std::cout << "   Epsilon    = " << fEpsilon   << std::endl;
  double  valsOut[6] = { 0., 0., 0., 0., 0., 0. };
  yTrackOutFT.DumpToArray(valsOut);
  endPosition=  ThreeVector( valsOut[0], valsOut[1], valsOut[2] );
  ThreeVector
  endMomentum=  ThreeVector( valsOut[3], valsOut[4], valsOut[3] );
  endDirection= invMomentumMag * endMomentum;
  // Check that endDirection is a unit vector here or later ?
#else
  //-------------------------------------------------------------------------------------
   // Do the single-track work using the Scalar Driver

  assert( fScalarDriver );
  if( verbose && !infoPrinted ) {
     std::cout << methodName << " > Using ScalarDriver " << fScalarDriver << std::endl;
     infoPrinted = true;
  }
  
  ScalarFieldTrack yTrackIn( startPosition, 
                             startDirection * startMomentumMag,
                             charge, 
                             0.0); // s_0  xo
  ScalarFieldTrack yTrackOut( yTrackIn );

  assert( fScalarDriver );
  // Simple call
  goodAdvance =
     fScalarDriver->AccurateAdvance( yTrackIn, step, fEpsilon, yTrackOut ); // , hInitial );
  endPosition=  yTrackOut.GetPosition();
  endDirection= yTrackOut.GetMomentumDirection();
#endif

  return goodAdvance;
}

VVectorField* GUFieldPropagator::GetField() 
{
   VVectorField* pField = nullptr;
   auto driver= GetScalarIntegrationDriver();
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

