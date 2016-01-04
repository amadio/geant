//
//
// Started from testPropagateMagField 85845 2014-11-05 15:43:58Z gcosmo $
//   Locate & Step within simple boxlike geometry, both
//   with and without voxels. Parameterised volumes are included.

//  UNFINISHED - adaptation from Geant4 test is incomplete
//   TODOs:
//     i) use either  VecGeom (or TGeo)
//    ii) use similar logic in stepping as GeantV (or proposed logic)
#include <assert.h>
// #include "ApproxEqual.hh"

// Global defs
// #include "globals.hh"
// #include "GUVPhysicalConstants.hh"
// #include "GUVSystemOfUnits.hh"

#include "ThreeVector.hh"
#include "base/Vector3D.h"
typedef vecgeom::Vector3D<double>  ThreeVector;

#include <iomanip>

// Build simple geometry:
// 4 small cubes + 1 slab (all GUVBoxes) are positioned inside a larger cuboid

#include "GUUniformMagField.hh"

#include "GUVExactHelixStepper.hh"
#include "GUVExplicitEuler.hh"
#include "GUVImplicitEuler.hh"
#include "GUVSimpleRunge.hh"
#include "GUVSimpleHeum.hh"
#include "GUVClassicalRK4.hh"
#include "GUVMag_UsualEqRhs.hh"
#include "GUVCashKarpRKF45.hh"
#include "GUVRKG3_Stepper.hh"
#include "GUVConstRK4.hh"
#include "GUVNystromRK4.hh"
#include "GUVHelixMixedStepper.hh"

#include "globals.hh"

GUVUniformMagField      uniformMagField(10.*tesla, 0., 0.); 
// GUVCachedMagneticField  myMagField( &uniformMagField, 1.0 * cm); 
// GUVString   fieldName("Uniform 10Tesla"); 

GUVQuadrupoleMagField   quadrupoleMagField( 10.*tesla/(50.*cm) ); 
GUVCachedMagneticField  myMagField( &quadrupoleMagField, 1.0 * cm); 
GUVString   fieldName("Cached Quadropole field, 20T/meter, cache=1cm"); 

GUVFieldManager* SetupField(int type)
{
    GUVFieldManager   *pFieldMgr;
    GUVChordFinder    *pChordFinder;
    GUVMag_UsualEqRhs *fEquation = new GUVMag_UsualEqRhs(&myMagField); 
    GUVMagIntegratorStepper *pStepper;

    std::cout << " Setting up field of type: " << fieldName << std::endl;

    switch ( type ) 
    {
      case 0: pStepper = new GUVExplicitEuler( fEquation ); break;
      case 1: pStepper = new GUVImplicitEuler( fEquation ); break;
      case 2: pStepper = new GUVSimpleRunge( fEquation ); break;
      case 3: pStepper = new GUVSimpleHeum( fEquation ); break;
      case 4: pStepper = new GUVClassicalRK4( fEquation ); break;
      case 5: pStepper = new GUVHelixExplicitEuler( fEquation ); break;
      case 6: pStepper = new GUVHelixImplicitEuler( fEquation ); break;
      case 7: pStepper = new GUVHelixSimpleRunge( fEquation ); break;
      case 8: pStepper = new GUVCashKarpRKF45( fEquation );    break;
      case 9: pStepper = new GUVExactHelixStepper( fEquation );   break;
      case 10: pStepper = new GUVRKG3_Stepper( fEquation );       break;
      case 11: pStepper = new GUVHelixMixedStepper( fEquation );  break;
      case 12: pStepper = new GUVConstRK4( fEquation ); break;
      case 13: pStepper = new GUVNystromRK4( fEquation ); break; 
      default: 
          pStepper = 0;   // Can use default= new GUVClassicalRK4( fEquation );
          GUVExceptionDescription ErrorMsg;
          ErrorMsg << " Incorrect Stepper type requested. Value was id= " 
                   << type << std::endl;
          ErrorMsg << " NO replacement stepper chosen! " << std::endl;
          GUVException("application::SetupField",
                      "Runtime Error",
                      FatalErrorInArgument,       //  use JustWarning,
                      " Invalid value of stepper type" );
          break; 
    }
    
    pFieldMgr= GUVTransportationManager::GetTransportationManager()->
       GetFieldManager();

    pFieldMgr->SetDetectorField( &myMagField );

    pChordFinder = new GUVChordFinder( &myMagField,
				      1.0e-2 * mm,
				      pStepper);
    pChordFinder->SetVerbose(1);  // ity();

    pFieldMgr->SetChordFinder( pChordFinder );

    return    pFieldMgr;
}

#include "GUVSimpleLocator.hh"
#include "GUVBrentLocator.hh"
#include "GUVMultiLevelLocator.hh"

GUVPropagatorInField*  SetupPropagator( int type)
{
    // GUVFieldManager* fieldMgr= 
    SetupField( type) ;

    // GUVChordFinder  theChordFinder( &MagField, 0.05*mm ); // Default stepper
 
    GUVPropagatorInField *thePropagator = 
      GUVTransportationManager::GetTransportationManager()->
       GetPropagatorInField ();

    // Let us test the new Minimum Epsilon Step functionality
    // thePropagator -> SetMinimumEpsilonStep( 1.0e-3 ) ; 
    // thePropagator -> SetMaximumEpsilonStep( 1.0e-5 ) ; 

    GUVNavigator *theNavigator= GUVTransportationManager::GetTransportationManager()->
       GetNavigatorForTracking();
    // Test the options for Locator
    GUVVIntersectionLocator *pLocator=0;
    std::cout << "Over-riding  PropagatorInField to use ";
    pLocator= new GUVMultiLevelLocator(theNavigator); std::cout << "Multi"; // default
    // pLocator= new GUVSimpleLocator(theNavigator); std::cout << "Simple";
    // pLocator= new GUVBrentLocator(theNavigator); std::cout << " Brent "; 
    std::cout << " Locator. ( In the unit test code. ) " << std::endl;

    thePropagator->SetIntersectionLocator(pLocator);

    return thePropagator;
}

GUVPropagatorInField *pMagFieldPropagator=0; 
//
// Test Stepping
//
bool testGUVPropagatorInField(GUVVPhysicalVolume*,     // *pTopNode, 
			       int             type)
{
    void report_endPV(GUVThreeVector    Position, 
                  GUVThreeVector UnitVelocity,
		  double step_len, 
                  double physStep, 
                  double safety,
		  GUVThreeVector EndPosition, 
                  GUVThreeVector EndUnitVelocity,
                  int             Step, 
                  GUVVPhysicalVolume* startVolume);

    GUVUniformMagField MagField(10.*tesla, 0., 0.);
    GUVNavigator   *pNavig= GUVTransportationManager::
                    GetTransportationManager()-> GetNavigatorForTracking();
    
    pMagFieldPropagator= SetupPropagator(type);

    double particleCharge= +1.0;  // in e+ units
    double spin=0.0;              // ignore the spin
    double magneticMoment= 0.0;   // ignore the magnetic moment

    GUVChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0); 

    GUVEquationOfMotion* equationOfMotion = 
        ( pMagFieldPropagator->GetChordFinder()->GetIntegrationDriver()->GetStepper())
        ->GetEquationOfMotion();
    
    equationOfMotion->SetChargeMomentumMass( chargeState, 
			            0.5 * proton_mass_c2, // Momentum in Mev/c
					 proton_mass_c2 );
    // pNavig->SetWorldVolume(pTopNode);

    GUVVPhysicalVolume *located;
    double step_len, physStep, safety;
    GUVThreeVector xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    GUVThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    
    // physStep=kInfinity;
    GUVThreeVector Position(0.,0.,0.); 
    GUVThreeVector UnitMomentum(0.,0.6,0.8);  
    GUVThreeVector EndPosition, EndUnitMomentum;

//
// Test location & Step computation
//  
    /* assert(located->GetName()=="World"); */
    if( std::fabs(UnitMomentum.mag() - 1.0) > 1.e-8 ) 
    {
      std::cerr << "UnitMomentum.mag() - 1.0 = " << UnitMomentum.mag() - 1.0 <<
	std::endl;
    }

    std::cout << std::endl; 

    for( int iparticle=0; iparticle < 2; iparticle++ )
    { 
       physStep=  2.5 * mm ;  // millimeters 
       Position = GUVThreeVector(0.,0.,0.) 
	        + iparticle * GUVThreeVector(0.2, 0.3, 0.4); 
       UnitMomentum = (GUVThreeVector(0.,0.6,0.8) 
		    + (float)iparticle * GUVThreeVector(0.1, 0.2, 0.3)).unit();

       double momentum = (0.5+iparticle*10.0) * proton_mass_c2; 

       double kineticEnergy =  momentum*momentum /
                  ( std::sqrt( momentum*momentum + proton_mass_c2 * proton_mass_c2 ) 
		    + proton_mass_c2 );
       double velocity = momentum / ( proton_mass_c2 + kineticEnergy );
       double labTof= 10.0*ns, properTof= 0.1*ns;
       GUVThreeVector Spin(1.0, 0.0, 0.0);
                                                   // Momentum in Mev/c ?

       GUVChargeState chargeSt(1.0, 0.0, 0.5 ); 
       // pMagFieldPropagator
       equationOfMotion->SetChargeMomentumMass(
		      chargeSt,                    // charge in e+ units
		      momentum, 
		      proton_mass_c2); 
       std::cout << std::endl;
       std::cout << "Test PropagateMagField: ***********************" << std::endl
            << " Starting New Particle with Position " << Position << std::endl 
	    << " and UnitVelocity " << UnitMomentum << std::endl;
       std::cout << " Momentum in GeV/c is " << momentum / GeV
	      << " = " << (0.5+iparticle*10.0)*proton_mass_c2 / MeV << " MeV"
              << std::endl;


       for( int istep=0; istep < 14; istep++ ){ 
          // std::cerr << "UnitMomentum Magnitude is " << UnitMomentum.mag() << std::endl;
	  located= pNavig->LocateGlobalPointAndSetup(Position);
	  // std::cerr << "Starting Step " << istep << " in volume " 
	       // << located->GetName() << std::endl;

          GUVFieldTrack  initTrack( Position, 
				   UnitMomentum,
				   0.0,            // starting S curve len
				   kineticEnergy,
				   proton_mass_c2,
				   velocity,
				   labTof, 
				   properTof,
				   0              // or &Spin
				   ); 

	  step_len=pMagFieldPropagator->ComputeStep( initTrack, 
						     physStep, 
						     safety,
						     located);
	  //       --------------------
	  EndPosition=     pMagFieldPropagator->EndPosition();
	  EndUnitMomentum= pMagFieldPropagator->EndMomentumDir();
	  //       --------
	  
	  if( std::fabs(EndUnitMomentum.mag2() - 1.0) > 1.e-8 )
	    std::cerr << "EndUnitMomentum.mag2() - 1.0 = " <<
	      EndUnitMomentum.mag2() - 1.0 << std::endl;

	  GUVThreeVector MoveVec = EndPosition - Position;
	  assert( MoveVec.mag() < physStep*(1.+1.e-9) );

	  // std::cout << " testPropagatorInField: After stepI " << istep  << " : " << std::endl;
	  report_endPV(Position, UnitMomentum, step_len, physStep, safety,
		       EndPosition, EndUnitMomentum, istep, located );

	  assert(safety>=0);
	  pNavig->SetGeometricallyLimitedStep();
	  // pMagFieldPropagator->SetGeometricallyLimitedStep();

	  Position= EndPosition;
	  UnitMomentum= EndUnitMomentum;
	  physStep *= 2.; 
       } // ...........................  end for ( istep )

       myMagField.ReportStatistics(); 

    }    // ..............................  end for ( iparticle )

    return(1);
}

void report_endPV(GUVThreeVector    Position, 
                  GUVThreeVector    InitialUnitVelocity,
		  double step_len, 
                  double physStep, 
                  double safety,
		  GUVThreeVector EndPosition, 
                  GUVThreeVector EndUnitVelocity,
                  int             Step, 
                  GUVVPhysicalVolume* startVolume)
		  //   GUVVPhysicalVolume* endVolume)
{
    const int verboseLevel=1;
    
    if( Step == 0 && verboseLevel <= 3 )
    {
       std::cout.precision(6);
       // std::cout.setf(ios_base::fixed,ios_base::floatfield);
       std::cout << std::setw( 5) << "Step#" << " "
            << std::setw( 9) << "X(mm)" << " "
            << std::setw( 9) << "Y(mm)" << " "  
            << std::setw( 9) << "Z(mm)" << " "
            << std::setw( 9) << " N_x " << " "
            << std::setw( 9) << " N_y " << " "
            << std::setw( 9) << " N_z " << " "
            << std::setw( 9) << " Delta|N|" << " "
            << std::setw( 9) << " Delta(N_z) " << " "
	   // << std::setw( 9) << "KinE(MeV)" << " "
	   // << std::setw( 9) << "dE(MeV)" << " "  
            << std::setw( 9) << "StepLen" << " "  
            << std::setw( 9) << "PhsStep" << " "  
            << std::setw( 9) << "Safety" << " "  
            << std::setw(18) << "NextVolume" << " "
            << std::endl;
    }
    //
    //
    if( verboseLevel > 3 )
    {
       std::cout << "End  Position is " << EndPosition << std::endl 
	    << " and UnitVelocity is " << EndUnitVelocity << std::endl;
       std::cout << "Step taken was " << step_len  
	    << " out of PhysicalStep= " <<  physStep << std::endl;
       std::cout << "Final safety is: " << safety << std::endl;

       std::cout << "Chord length = " << (EndPosition-Position).mag() << std::endl;
       std::cout << std::endl; 
    }
    else // if( verboseLevel > 0 )
    {
       std::cout.precision(6);
       std::cout << std::setw( 5) << Step << " "
	    << std::setw( 9) << Position.x() << " "
	    << std::setw( 9) << Position.y() << " "
	    << std::setw( 9) << Position.z() << " "
	    << std::setw( 9) << EndUnitVelocity.x() << " "
	    << std::setw( 9) << EndUnitVelocity.y() << " "
	      << std::setw( 9) << EndUnitVelocity.z() << " ";
       std::cout.precision(2); 
       std::cout
	    << std::setw( 9) << EndUnitVelocity.mag()-InitialUnitVelocity.mag() << " "
	    << std::setw( 9) << EndUnitVelocity.z() - InitialUnitVelocity.z() << " ";
	 //    << std::setw( 9) << KineticEnergy << " "
	 //    << std::setw( 9) << EnergyDifference << " "
       std::cout.precision(6);
       std::cout 
	    << std::setw( 9) << step_len << " "
	    << std::setw( 9) << physStep << " "
	    << std::setw( 9) << safety << " ";
       if( startVolume != 0) {
	 std::cout << std::setw(12) << startVolume->GetName() << " ";
       } else {
	 std::cout << std::setw(12) << "OutOfWorld" << " ";
       }
#if 0
       if( endVolume != 0) 
	 std::cout << std::setw(12) << endVolume()->GetName() << " ";
       else 
	 std::cout << std::setw(12) << "OutOfWorld" << " ";
#endif
       std::cout << std::endl;
    }
}

// Main program
// -------------------------------
int main(int argc, char **argv)
{
    GUVVPhysicalVolume *myTopNode;
    int type, optim, optimSaf;
    bool optimiseVoxels=true;
    bool optimisePiFwithSafety=true;

    type = 8 ;
    std::cout << " Arguments:  stepper-no  optimise-Voxels optimise-PiF-with-safety" << std::endl;

    if( argc >= 2 ){
       type = atoi(argv[1]);
    } 

    if( argc >=3 ){
      optim= atoi(argv[2]);
      if( optim == 0 ) { optimiseVoxels = false; }
    }

    if( argc >=4 ){
      optimSaf= atoi(argv[3]);
      if( optimSaf == 0 ) { optimisePiFwithSafety= false; }
    }

    std::cout << " Testing with stepper number    " << type << std::endl; 
    std::cout << "             " ; 
    std::cout << " voxel optimisation      " ; 
    // if (optimiseVoxels)   std::cout << "On"; 
    // else                  std::cout << "Off"; 
    std::cout << (optimiseVoxels ? "On" : "Off")  << std::endl;
    std::cout << "             " ; 
    std::cout << " Propagator safety optim " ; 
    // const char* OnOff= (optimisePiFwithSafety ? "on" : "off") ; 
    // std::cout << OnOff << std::endl;
    std::cout << (optimisePiFwithSafety ? "On" : "Off")  << std::endl;

    // Create the geometry & field 
    myTopNode=BuildGeometry();	// Build the geometry
 
    GUVNavigator *pNavig= GUVTransportationManager::
                    GetTransportationManager()-> GetNavigatorForTracking();
    pNavig->SetWorldVolume(myTopNode);

    GUVGeometryManager::GetInstance()->CloseGeometry(false);

    // Setup the propagator (will be overwritten by testGUVPropagator ...)
    pMagFieldPropagator= SetupPropagator(type);
    std::cout << " Using default values for " 
	   << " Min Eps = "  <<   pMagFieldPropagator->GetMinimumEpsilonStep()
           << " and "
	   << " MaxEps = " <<  pMagFieldPropagator->GetMaximumEpsilonStep()
	   << std::endl; 

    pMagFieldPropagator->SetUseSafetyForOptimization(optimisePiFwithSafety); 

// Do the tests without voxels
    std::cout << " Test with no voxels" << std::endl; 
    testGUVPropagatorInField(myTopNode, type);

    pMagFieldPropagator->SetUseSafetyForOptimization(optimiseVoxels); 
    pMagFieldPropagator->SetVerboseLevel( 1 ); 

// Repeat tests but with full voxels
    std::cout << " Test with full voxels" << std::endl; 

    GUVGeometryManager::GetInstance()->OpenGeometry();
    GUVGeometryManager::GetInstance()->CloseGeometry(true);

    testGUVPropagatorInField(myTopNode, type);

    GUVGeometryManager::GetInstance()->OpenGeometry();

    std::cout << std::endl
	   << "----------------------------------------------------------"
	   << std::endl; 

// Repeat tests with full voxels and modified parameters
    std::cout << "Test with more accurate parameters " << std::endl; 

    double  maxEpsStep= 0.001;
    double  minEpsStep= 2.5e-8;
    std::cout << " Setting values for Min Eps = " << minEpsStep 
           << " and MaxEps = " << maxEpsStep << std::endl; 

    pMagFieldPropagator->SetMaximumEpsilonStep(maxEpsStep);
    pMagFieldPropagator->SetMinimumEpsilonStep(minEpsStep);

    GUVGeometryManager::GetInstance()->OpenGeometry();
    GUVGeometryManager::GetInstance()->CloseGeometry(true);

    testGUVPropagatorInField(myTopNode, type);

    GUVGeometryManager::GetInstance()->OpenGeometry();

    optimiseVoxels = ! optimiseVoxels;
// Repeat tests but with the opposite optimisation choice
    std::cout << " Now test with optimisation " ; 
    if (optimiseVoxels)   std::cout << "on"; 
    else            std::cout << "off"; 
    std::cout << std::endl;

    pMagFieldPropagator->SetUseSafetyForOptimization(optimiseVoxels); 
    testGUVPropagatorInField(myTopNode, type);

    GUVGeometryManager::GetInstance()->OpenGeometry();

    // Cannot delete GUVTransportationManager::GetInstance(); 
    return 0;
}


  
