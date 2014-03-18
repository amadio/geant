//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// 
// 
//  This class implements an algorithm to track a particle in a
//  non-uniform magnetic field. It utilises an ODE solver (with
//  the Runge - Kutta method) to evolve the particle, and drives it
//  until the particle has traveled a set distance or it enters a new 
//  volume.
//                                                                     
// 14.10.96 John Apostolakis,   design and implementation
// 17.03.97 John Apostolakis,   renaming new set functions being added
//
// $Id: G4PropagatorInField.cc,v 1.52 2010-07-13 15:59:42 gcosmo Exp $
// GEANT4 tag $ Name:  $
// ---------------------------------------------------------------------------

//#include "G4ios.hh"
//#include <iomanip>
//#include "G4Navigator.hh"
//#include "G4GeometryTolerance.hh"
//#include "G4VPhysicalVolume.hh"
//#include "G4MultiLevelLocator.hh"
//#include "G4VCurvedTrajectoryFilter.hh"

#include "GPPropagatorInField.h"
#include "GPChordFinder.h"
#include "GPConstants.h"
#include "GPUtils.h"

///////////////////////////////////////////////////////////////////////////
//
// Constructors and destructor

FQUALIFIER
void GPPropagatorInField_Constructor( GPPropagatorInField *This, 
				      GPNavigator    *theNavigator, 
				      GPFieldManager *detectorFieldMgr,
				      GPMultiLevelLocator *vLocator )
{
  This->fMax_loop_count = 1000;
  This->fUseSafetyForOptimisation = true; 
  //(false) is less sensitive to incorrect safety
  This->fZeroStepThreshold = 0.0; //length of what is recognised as 'zero' step

  This->fDetectorFieldMgr = detectorFieldMgr;
  //  fpTrajectoryFilter = 0;

  This->fNavigator = theNavigator;

  This->fCurrentFieldMgr = detectorFieldMgr;

  This->fSetFieldMgr = false;
  This->fCharge = 0.0;
  This->fInitialMomentumModulus = 0.0;
  This->fMass = 0.0;

  GPFieldTrack_Constructor2(&(This->End_PointAndTangent),
			    GPThreeVector_create(0.,0.,0.),
			    GPThreeVector_create(0.,0.,0.),
			    0.0,0.0,0.0,0.0,0.0,
			    0);

  This->fParticleIsLooping = false;
  This->fNoZeroStep = 0;
  This->fVerboseLevel = 0;

  if(This->fDetectorFieldMgr) 
  { 
    This->fEpsilonStep = 
      GPFieldManager_GetMaximumEpsilonStep(This->fDetectorFieldMgr);
  }
  else                  
  { 
    This->fEpsilonStep= 1.0e-5; 
  } 

  This->fActionThreshold_NoZeroSteps = 2; 
  This->fSevereActionThreshold_NoZeroSteps = 10; 
  This->fAbandonThreshold_NoZeroSteps = 50; 
  This->fFull_CurveLen_of_LastAttempt = -1; 
  This->fLast_ProposedStepLength = -1;
  This->fLargestAcceptableStep = 1000.0 * meter;

  This->fPreviousSftOrigin = GPThreeVector_create(0.,0.,0.);
  This->fPreviousSafety= 0.0;

  //  This->kCarTolerance = 
  //    G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  This->kCarTolerance = 1E-9*millimeter;

  This->fZeroStepThreshold= 
    GPfmax( 1.0e5 * This->kCarTolerance, 1.0e-1 * micrometer ) ; 

  // Definding Intersection Locator and his parameters

  if(vLocator==0){
    //    This->fIntersectionLocator= new G4MultiLevelLocator(theNavigator);
    GPMultiLevelLocator mLocator;
    GPMultiLevelLocator_Constructor(&mLocator,theNavigator);
    This->fIntersectionLocator = &mLocator;
    This->fAllocatedLocator=true;
  }else{
    This->fIntersectionLocator=vLocator;
    This->fAllocatedLocator=false;
  }

  GPPropagatorInField_RefreshIntersectionLocator(This);//Copy all relevant parameters 
}

//GPPropagatorInField::~GPPropagatorInField()
//{
//  if(fAllocatedLocator)delete  fIntersectionLocator; 
//}

//********************************************************************
// Update the IntersectionLocator with current parameters

FQUALIFIER
void GPPropagatorInField_RefreshIntersectionLocator(GPPropagatorInField *This)
{
  GPVIntersectionLocator_SetEpsilonStepFor(This->fIntersectionLocator,
					   This->fEpsilonStep);
  GPVIntersectionLocator_SetDeltaIntersectionFor(This->fIntersectionLocator,
		 GPFieldManager_GetDeltaIntersection(This->fCurrentFieldMgr));
  GPVIntersectionLocator_SetChordFinderFor(This->fIntersectionLocator,
			 GPPropagatorInField_GetChordFinder(This));
  GPVIntersectionLocator_SetSafetyParametersFor(This->fIntersectionLocator, 
					       This->fUseSafetyForOptimisation);
}

///////////////////////////////////////////////////////////////////////////
//
// Compute the next geometric Step

FQUALIFIER
G4double
GPPropagatorInField_ComputeStep(GPPropagatorInField *This,
				GPFieldTrack&      pFieldTrack,
				G4double           CurrentProposedStepLength,
				G4double&          currentSafety,   // IN/OUT
				GPVPhysicalVolume* pPhysVol ) 
{  

  // If CurrentProposedStepLength is too small for finding Chords
  // then return with no action (for now - TODO: some action)
  //
  if(CurrentProposedStepLength < This->kCarTolerance)
  {
    return kInfinity;
  }

  // Introducing smooth trajectory display (jacek 01/11/2002)
  //
  //  if (fpTrajectoryFilter)
  //  {
  //    fpTrajectoryFilter->CreateNewTrajectorySegment();
  //  }

  // Parameters for adaptive Runge-Kutta integration
  
  G4double      h_TrialStepSize;        // 1st Step Size 
  G4double      TruePathLength = CurrentProposedStepLength;
  G4double      StepTaken = 0.0; 
  G4double      s_length_taken, epsilon ; 
  G4bool        intersects;
  G4bool        first_substep = true;

  G4double      NewSafety;
  This->fParticleIsLooping = false;

  // If not yet done, 
  //   Set the field manager to the local  one if the volume has one, 
  //                      or to the global one if not
  //
  //*************************************************************************
  //@@@G4FWP the field map is directly used from a bfield map
  //  if( !fSetFieldMgr ) fCurrentFieldMgr= FindAndSetFieldManager( pPhysVol ); 
  if( !(This->fSetFieldMgr) ) 
    This->fCurrentFieldMgr = GPPropagatorInField_FindAndSetFieldManager(This); 
  //@@@G4FWP*****************************************************************
  // For the next call, the field manager must again be set
  This->fSetFieldMgr= false;

  // GetChordFinder()->SetChargeMomentumMass(fCharge, 
  //                                         fInitialMomentumModulus, fMass); 
  GPChordFinder *chordFinder = 
    GPFieldManager_GetChordFinder(This->fCurrentFieldMgr);
  GPChordFinder_SetChargeMomentumMass(chordFinder,
				      This->fCharge, 
				      This->fInitialMomentumModulus, 
				      This->fMass);

  // Values for Intersection Locator has to be updated on each call for the
  // case that CurrentFieldManager has changed from the one of previous step
  GPPropagatorInField_RefreshIntersectionLocator(This);

  //  G4FieldTrack  CurrentState(pFieldTrack);
  GPFieldTrack&  CurrentState  = pFieldTrack;
  //  G4FieldTrack  OriginalState = CurrentState;
  GPFieldTrack&  OriginalState = pFieldTrack;

  // If the Step length is "infinite", then an approximate-maximum Step
  // length (used to calculate the relative accuracy) must be guessed.
  //
  if( CurrentProposedStepLength >= This->fLargestAcceptableStep )
  {
    GPThreeVector StartPointA, VelocityUnit;
    StartPointA  = GPFieldTrack_GetPosition(&pFieldTrack);
    VelocityUnit = GPFieldTrack_GetMomentumDir(&pFieldTrack);

    G4double trialProposedStep = 1.e2 * ( 10.0 * centimeter + 
    //   	fNavigator->GetWorldVolume()->GetLogicalVolume()->
    //          GetSolid()->DistanceToOut(StartPointA, VelocityUnit) );
      GPVSolid_DistanceToOut2(
	GPLogicalVolume_GetSolid(
	  GPVPhysicalVolume_GetLogicalVolume(
	    GPNavigator_GetWorldVolume(This->fNavigator))),
	StartPointA, VelocityUnit,0,0,0));

    CurrentProposedStepLength= GPfmin( trialProposedStep,
				       This->fLargestAcceptableStep ); 
  }
  epsilon = 
    GPFieldManager_GetDeltaOneStep(This->fCurrentFieldMgr)/CurrentProposedStepLength;

  // G4double raw_epsilon= epsilon;
  G4double epsilonMin= GPFieldManager_GetMinimumEpsilonStep(This->fCurrentFieldMgr);
  G4double epsilonMax= GPFieldManager_GetMaximumEpsilonStep(This->fCurrentFieldMgr);
  if( epsilon < epsilonMin ) epsilon = epsilonMin;
  if( epsilon > epsilonMax ) epsilon = epsilonMax;
  GPPropagatorInField_SetEpsilonStep(This, epsilon );

  // G4cout << "G4PiF: Epsilon of current step - raw= " << raw_epsilon
  //        << " final= " << epsilon << G4endl;

  // Shorten the proposed step in case of earlier problems (zero steps)
  // 

  if( This->fNoZeroStep > This->fActionThreshold_NoZeroSteps )
  {
    G4double stepTrial;

    stepTrial= This->fFull_CurveLen_of_LastAttempt; 
    if( (stepTrial <= 0.0) && (This->fLast_ProposedStepLength > 0.0) ) 
      stepTrial= This->fLast_ProposedStepLength; 

    G4double decreaseFactor = 0.9; // Unused default
    if(   (This->fNoZeroStep < This->fSevereActionThreshold_NoZeroSteps)
       && (stepTrial > 100.0*This->fZeroStepThreshold) )
    {
      // Attempt quick convergence
      //
      decreaseFactor= 0.25;
    } 
    else
    {
      // We are in significant difficulties, probably at a boundary that
      // is either geometrically sharp or between very different materials.
      // Careful decreases to cope with tolerance are required.
      //
      if( stepTrial > 100.0*This->fZeroStepThreshold )
        decreaseFactor = 0.35;     // Try decreasing slower
      else if( stepTrial > 100.0*This->fZeroStepThreshold )
        decreaseFactor= 0.5;       // Try yet slower decreases
      else if( stepTrial > 10.0*This->fZeroStepThreshold )
        decreaseFactor= 0.75;      // Try even slower decreases
      else
        decreaseFactor= 0.9;       // Try very slow decreases
     }
     stepTrial *= decreaseFactor;

     if( stepTrial == 0.0 )  //  Change to make it < 0.1 * kCarTolerance ??
     {
       //     std::ostringstream message;
       //     message << "Particle abandoned due to lack of progress in field."
       //               << G4endl
       //        << "  Properties : " << pFieldTrack << G4endl
       //        << "  Attempting a zero step = " << stepTrial << G4endl
       //        << "  while attempting to progress after " << fNoZeroStep
       //        << " trial steps. Will abandon step.";
       //G4Exception("GPPropagatorInField_ComputeStep()", "GeomNav1002",
       //            JustWarning, message);
       This->fParticleIsLooping= true;
       return 0;  // = stepTrial;
     }
     if( stepTrial < CurrentProposedStepLength )
       CurrentProposedStepLength = stepTrial;
  }
  This->fLast_ProposedStepLength = CurrentProposedStepLength;

  G4int do_loop_count = 0; 
  do
  { 
    GPFieldTrack& SubStepStartState = CurrentState;

    GPThreeVector SubStartPoint = GPFieldTrack_GetPosition(&CurrentState); 

    if( !first_substep) {
      GPNavigator_LocateGlobalPointWithinVolume(This->fNavigator,SubStartPoint);
    }

    // How far to attempt to move the particle !
    //
    h_TrialStepSize = CurrentProposedStepLength - StepTaken;

    // Integrate as far as "chord miss" rule allows.
    //

    s_length_taken = GPChordFinder_AdvanceChordLimited(
                             GPPropagatorInField_GetChordFinder(This), 
                             CurrentState,    // Position & velocity
                             h_TrialStepSize,
                             This->fEpsilonStep,
                             This->fPreviousSftOrigin,
                             This->fPreviousSafety
                             );
    //  CurrentState is now updated with the final position and velocity. 

    This->fFull_CurveLen_of_LastAttempt = s_length_taken;

    GPThreeVector  EndPointB = GPFieldTrack_GetPosition(&CurrentState); 
    GPThreeVector  InterSectionPointE;
    G4double       LinearStepLength;
 
    // Intersect chord AB with geometry
    intersects= GPPropagatorInField_IntersectChord(This, 
    						   SubStartPoint, EndPointB, 
    						   NewSafety, LinearStepLength, 
    						   InterSectionPointE );
    // E <- Intersection Point of chord AB and either volume A's surface 
    //                                  or a daughter volume's surface ..

    NewSafety = 0;
    if( first_substep ) { 
       currentSafety = NewSafety;
    } // Updating safety in other steps is potential future extention

    if( intersects )
    {
      GPFieldTrack& IntersectPointVelct_G = CurrentState;  // FT-Def-Construct

      // Find the intersection point of AB true path with the surface
      //   of vol(A), if it exists. Start with point E as first "estimate".
      G4bool recalculatedEndPt= false;

      G4bool found_intersection = 
	GPMultiLevelLocator_EstimateIntersectionPoint(
                                    This->fIntersectionLocator,
				    SubStepStartState, CurrentState, 
				    InterSectionPointE, IntersectPointVelct_G,
				    recalculatedEndPt,This->fPreviousSafety,
				    This->fPreviousSftOrigin);

      intersects = intersects && found_intersection;
      if( found_intersection ) {        
	This->End_PointAndTangent= IntersectPointVelct_G;  //G is our EndPoint
	StepTaken = TruePathLength = 
	    GPFieldTrack_GetCurveLength(&IntersectPointVelct_G)
	  - GPFieldTrack_GetCurveLength(&OriginalState);
      } else {
	// intersects= false;          // "Minor" chords do not intersect
	if( recalculatedEndPt ){
	  CurrentState= IntersectPointVelct_G; 
	}
      }
    }

    if( !intersects )
    {
      StepTaken += s_length_taken; 
      // For smooth trajectory display (jacek 01/11/2002)
      //      if (This->fpTrajectoryFilter) {
      // fpTrajectoryFilter->TakeIntermediatePoint(CurrentState.GetPosition());
      //      }
    }
    first_substep = false;

    do_loop_count++;

  } while( (!intersects )
        && (StepTaken + This->kCarTolerance < CurrentProposedStepLength)  
        && ( do_loop_count < This->fMax_loop_count ) );

  if( do_loop_count >= This->fMax_loop_count  )
  {
    This->fParticleIsLooping = true;
  }

  if( !intersects )
  {
    // Chord AB or "minor chords" do not intersect
    // B is the endpoint Step of the current Step.
    //
    This->End_PointAndTangent = CurrentState; 
    TruePathLength = StepTaken;
  }
  
  // Set pFieldTrack to the return value
  //
  pFieldTrack = This->End_PointAndTangent;

  // In particular anomalous cases, we can get repeated zero steps
  // In order to correct this efficiently, we identify these cases
  // and only take corrective action when they occur.
  // 
  if( ( (TruePathLength < This->fZeroStepThreshold) 
	&& ( TruePathLength+This->kCarTolerance < CurrentProposedStepLength  ) 
	) 
      || ( TruePathLength < 0.5*This->kCarTolerance )
    )
  {
    This->fNoZeroStep++;
  }
  else{
    This->fNoZeroStep = 0;
  }

  if( This->fNoZeroStep > This->fAbandonThreshold_NoZeroSteps )
  { 
     This->fParticleIsLooping = true;
     //     std::ostringstream message;
     //     message << "Particle is stuck; it will be killed." << G4endl
     //     << "  Zero progress for "  << fNoZeroStep << " attempted steps." 
     //     << G4endl
     //     << "  Proposed Step is " << CurrentProposedStepLength
     //     << " but Step Taken is "<< fFull_CurveLen_of_LastAttempt << G4endl
     //     << "  For Particle with Charge = " << fCharge
     //     << " Momentum = "<< fInitialMomentumModulus
     //     << " Mass = " << fMass << G4endl;
     //if( pPhysVol )
     //  message << " in volume " << pPhysVol->GetName() ; 
     //else
     //  message << " in unknown or null volume. " ; 
     //G4Exception("GPPropagatorInField_ComputeStep()",
     //            "GeomNav1002", JustWarning, message);
     This->fNoZeroStep = 0; 
  }
 
  return TruePathLength;
}

///////////////////////////////////////////////////////////////////////////
//
// Dumps status of propagator.

FQUALIFIER
void
GPPropagatorInField_printStatus( const GPFieldTrack&,   //     StartFT,
				 const GPFieldTrack&,   //     CurrentFT, 
				 G4double,             //requestStep, 
				 G4double,             //safety,
				 G4int                //stepNo
				 // G4VPhysicalVolume*   startVolume
				 )
{
  //do nothing for GPU
  ;
}

///////////////////////////////////////////////////////////////////////////
//
// Prints Step diagnostics

FQUALIFIER
void GPPropagatorInField_PrintStepLengthDiagnostic(
			 G4double CurrentProposedStepLength,
			 G4double decreaseFactor,
			 G4double stepTrial,
			 const GPFieldTrack& )
{
  //do nothing for GPU
  ;
}

// Access the points which have passed through the filter. The
// points are stored as ThreeVectors for the initial impelmentation
// only (jacek 30/10/2002)
// Responsibility for deleting the points lies with
// SmoothTrajectoryPoint, which is the points' final
// destination. The points pointer is set to NULL, to ensure that
// the points are not re-used in subsequent steps, therefore THIS
// METHOD MUST BE CALLED EXACTLY ONCE PER STEP. (jacek 08/11/2002)

//std::vector<G4ThreeVector>*
FQUALIFIER
void GPPropagatorInField_GimmeTrajectoryVectorAndForgetIt() 
{
  //do nothing for GPU
  ;

  // NB, GimmeThePointsAndForgetThem really forgets them, so it can
  // only be called (exactly) once for each step.

  //  if (fpTrajectoryFilter)
  //  {
  //    return fpTrajectoryFilter->GimmeThePointsAndForgetThem();
  //  }
  //  else
  //  {
  //    return 0;
  //  }
}

//void 
//GPPropagatorInField_SetTrajectoryFilter(G4VCurvedTrajectoryFilter* filter)
//{
//  This->fpTrajectoryFilter = filter;
//}
//

FQUALIFIER
void GPPropagatorInField_ClearPropagatorState(GPPropagatorInField *This)
{
  // Goal: Clear all memory of previous steps,  cached information

  This->fParticleIsLooping= false;
  This->fNoZeroStep= 0;

  GPFieldTrack_Constructor2(&(This->End_PointAndTangent), 
			    GPThreeVector_create(0.,0.,0.),
			    GPThreeVector_create(0.,0.,0.),
			    0.0,0.0,0.0,0.0,0.0,
			    0); 
  This->fFull_CurveLen_of_LastAttempt = -1; 
  This->fLast_ProposedStepLength = -1;

  This->fPreviousSftOrigin= GPThreeVector_create(0.,0.,0.);
  This->fPreviousSafety= 0.0;
}

//G4FieldManager* 
//GPPropagatorInField_FindAndSetFieldManager( G4VPhysicalVolume* pCurrentPhysicalVolume)

FQUALIFIER
GPFieldManager* 
GPPropagatorInField_FindAndSetFieldManager(GPPropagatorInField *This)
{
  GPFieldManager* currentFieldMgr;

  currentFieldMgr = This->fDetectorFieldMgr;

  //@@@G4FWP skips this part since the magnetic field map is provided in a 
  // different way in GPU - so no need to know the current physical volume

  //  if( pCurrentPhysicalVolume)
  //  {
  //     G4FieldManager *pRegionFieldMgr= 0, *localFieldMgr = 0;
  //     G4LogicalVolume* pLogicalVol= pCurrentPhysicalVolume->GetLogicalVolume();
  //
  //     if( pLogicalVol ) { 
	// Value for Region, if any, Overrides 
  //	G4Region*  pRegion= pLogicalVol->GetRegion();
  //	if( pRegion ) { 
  //  //	   pRegionFieldMgr= pRegion->GetFieldManager();
  //	   if( pRegionFieldMgr ) 
  //	     currentFieldMgr= pRegionFieldMgr;
  //	}
  //
	// 'Local' Value from logical volume, if any, Overrides 
  //	localFieldMgr= pLogicalVol->GetFieldManager();
  //	if ( localFieldMgr ) 
  //	   currentFieldMgr = localFieldMgr;
  //    }
  //  }
  This->fCurrentFieldMgr = currentFieldMgr;

  // Flag that field manager has been set.
  This->fSetFieldMgr= true;

  return currentFieldMgr;
}

FQUALIFIER
G4int GPPropagatorInField_SetVerboseLevel(GPPropagatorInField *This, G4int level )
{
  G4int oldval= This->fVerboseLevel;
  This->fVerboseLevel= level;

  // Forward the verbose level 'reduced' to ChordFinder,
  // MagIntegratorDriver ... ? 
  //
  //  G4MagInt_Driver* integrDriver= GetChordFinder()->GetIntegrationDriver(); 
  //  integrDriver->SetVerboseLevel( fVerboseLevel - 2 );
  //  G4cout << "Set Driver verbosity to " << fVerboseLevel - 2 << G4endl;

  return oldval;
}


// ------------------------------------------------------------------------
// G4PropagatorInField.icc inline implementation

FQUALIFIER
GPChordFinder* GPPropagatorInField_GetChordFinder( GPPropagatorInField *This )
{
  // The "Chord Finder" of the current Field Mgr is used
  //    -- this could be of the global field manager
  //        or that of another, from the current volume 
  return GPFieldManager_GetChordFinder( This->fCurrentFieldMgr ); 
}

FQUALIFIER
void GPPropagatorInField_SetChargeMomentumMass( GPPropagatorInField *This,  
						G4double Charge,    // in e+ units
						G4double Momentum,  // in GeV/c 
						G4double Mass)      // in ? units
{
  // GetChordFinder()->SetChargeMomentumMass(Charge, Momentum, Mass);
  //  --> Not needed anymore, as it is done in ComputeStep for the 
  //       ChordFinder of the current step (which is known only then).
  This->fCharge = Charge;
  This->fInitialMomentumModulus = Momentum;
  This->fMass = Mass; 
}

//  Obtain the final space-point and velocity (normal) at the end of the Step
//

FQUALIFIER
GPThreeVector  GPPropagatorInField_EndPosition( GPPropagatorInField *This )
{
  return  GPFieldTrack_GetPosition(&(This->End_PointAndTangent)); 
}

FQUALIFIER
GPThreeVector  GPPropagatorInField_EndMomentumDir( GPPropagatorInField *This )
{
  return  GPFieldTrack_GetMomentumDir(&(This->End_PointAndTangent)); 
}

FQUALIFIER
G4double GPPropagatorInField_GetEpsilonStep( GPPropagatorInField *This )
{ 
  return This->fEpsilonStep; 
}

FQUALIFIER
void     GPPropagatorInField_SetEpsilonStep( GPPropagatorInField *This,
					     G4double newEps )
{
  This->fEpsilonStep=newEps;
}

FQUALIFIER
G4bool   GPPropagatorInField_IsParticleLooping( GPPropagatorInField *This )
{
  return This->fParticleIsLooping;
}

FQUALIFIER
G4int    GPPropagatorInField_GetMaxLoopCount( GPPropagatorInField *This )
{
  return This->fMax_loop_count;
}

FQUALIFIER
void GPPropagatorInField_SetMaxLoopCount( GPPropagatorInField *This,
					      G4int new_max ) 
{
  This->fMax_loop_count = new_max;
}

// #if 0
FQUALIFIER
G4double GPPropagatorInField_GetDeltaIntersection( GPPropagatorInField *This )
{
  return GPFieldManager_GetDeltaIntersection(This->fCurrentFieldMgr);
} 

FQUALIFIER
G4double GPPropagatorInField_GetDeltaOneStep( GPPropagatorInField *This )
{
  return GPFieldManager_GetDeltaOneStep(This->fCurrentFieldMgr);
}
// #endif 

FQUALIFIER
G4int GPPropagatorInField_GetVerboseLevel( GPPropagatorInField *This )
{
  return This->fVerboseLevel;
}
FQUALIFIER
G4int GPPropagatorInField_Verbose( GPPropagatorInField *This )
{
  return GPPropagatorInField_GetVerboseLevel(This);
}

FQUALIFIER
GPFieldTrack GPPropagatorInField_GetEndState( GPPropagatorInField *This )
{
  return This->End_PointAndTangent;
}

// Minimum for Relative accuracy of a Step in volumes of global field
FQUALIFIER 
G4double GPPropagatorInField_GetMinimumEpsilonStep( GPPropagatorInField *This )
{
  return GPFieldManager_GetMinimumEpsilonStep(This->fDetectorFieldMgr);
}

FQUALIFIER 
void      GPPropagatorInField_SetMinimumEpsilonStep(GPPropagatorInField *This,
						    G4double newEpsMin )
{
  GPFieldManager_SetMinimumEpsilonStep(This->fDetectorFieldMgr, newEpsMin);
}

// Maximum for Relative accuracy of any Step 
FQUALIFIER 
G4double  GPPropagatorInField_GetMaximumEpsilonStep( GPPropagatorInField *This )
{
  return GPFieldManager_GetMaximumEpsilonStep(This->fDetectorFieldMgr);
}

FQUALIFIER 
void      GPPropagatorInField_SetMaximumEpsilonStep(GPPropagatorInField *This,
						    G4double newEpsMax )
{
  GPFieldManager_SetMaximumEpsilonStep(This->fDetectorFieldMgr, newEpsMax );
}

FQUALIFIER
void GPPropagatorInField_SetLargestAcceptableStep(GPPropagatorInField *This,
						  G4double newBigDist )
{
  if( This->fLargestAcceptableStep>0.0 )
  {
    This->fLargestAcceptableStep = newBigDist;
  }
}

FQUALIFIER
G4double  GPPropagatorInField_GetLargestAcceptableStep( GPPropagatorInField *This )
{
  return This->fLargestAcceptableStep;
}

FQUALIFIER
GPFieldManager*  GPPropagatorInField_GetCurrentFieldManager( GPPropagatorInField *This )
{
  return This->fCurrentFieldMgr;
} 
FQUALIFIER
void GPPropagatorInField_SetThresholdNoZeroStep( GPPropagatorInField *This,
						 G4int noAct,
						 G4int noHarsh,
						 G4int noAbandon )
{
  if( noAct>0 ) 
    This->fActionThreshold_NoZeroSteps = noAct; 

  if( noHarsh > This->fActionThreshold_NoZeroSteps )
    This->fSevereActionThreshold_NoZeroSteps = noHarsh; 
  else
    This->fSevereActionThreshold_NoZeroSteps = 2*(This->fActionThreshold_NoZeroSteps+1);

  if( noAbandon > This->fSevereActionThreshold_NoZeroSteps+5 )
    This->fAbandonThreshold_NoZeroSteps = noAbandon; 
  else
    This->fAbandonThreshold_NoZeroSteps = 2*(This->fSevereActionThreshold_NoZeroSteps+3); 
}

FQUALIFIER
G4int GPPropagatorInField_GetThresholdNoZeroSteps( GPPropagatorInField *This,
						   G4int i )
{
   G4int t=0;
   if( i==0 )     { t = 3; }     // No of parameters
   else if (i==1) { t = This->fActionThreshold_NoZeroSteps; }
   else if (i==2) { t = This->fSevereActionThreshold_NoZeroSteps; }
   else if (i==3) { t = This->fAbandonThreshold_NoZeroSteps; }

   return t;
}

FQUALIFIER 
G4double  GPPropagatorInField_GetZeroStepThreshold( GPPropagatorInField *This )
{ 
  return This->fZeroStepThreshold; 
}

FQUALIFIER 
void      GPPropagatorInField_SetZeroStepThreshold(GPPropagatorInField *This,
						   G4double newLength )
{ 
  This->fZeroStepThreshold= newLength;
}

FQUALIFIER
void GPPropagatorInField_SetDetectorFieldManager(GPPropagatorInField *This,
				      GPFieldManager* newDetectorFieldManager)
{
  This->fDetectorFieldMgr = newDetectorFieldManager; 
}

FQUALIFIER
void GPPropagatorInField_SetUseSafetyForOptimization(GPPropagatorInField *This,
						     G4bool value )
{
  This->fUseSafetyForOptimisation= value;
}

FQUALIFIER 
G4bool 
GPPropagatorInField_GetUseSafetyForOptimization( GPPropagatorInField *This ) 
{ 
  return This->fUseSafetyForOptimisation; 
}

FQUALIFIER 
void GPPropagatorInField_SetNavigatorForPropagating(GPPropagatorInField *This, 
					 GPNavigator *SimpleOrMultiNavigator )
{
  if(SimpleOrMultiNavigator)  { 
     This->fNavigator= SimpleOrMultiNavigator; 
     if( This->fIntersectionLocator ) {
       GPVIntersectionLocator_SetNavigatorFor(This->fIntersectionLocator, 
					      SimpleOrMultiNavigator );
     }
  }
}

FQUALIFIER
GPNavigator* GPPropagatorInField_GetNavigatorForPropagating( GPPropagatorInField *This )
{
  return This->fNavigator;
} 

FQUALIFIER 
void GPPropagatorInField_SetIntersectionLocator(GPPropagatorInField *This,
						GPMultiLevelLocator *pIntLoc )
{
  if(pIntLoc)  { This->fIntersectionLocator= pIntLoc; }
}

//G4VIntersectionLocator* GPPropagatorInField_GetIntersectionLocator( GPPropagatorInField *This )
FQUALIFIER
GPMultiLevelLocator* GPPropagatorInField_GetIntersectionLocator( GPPropagatorInField *This )
{
  return This->fIntersectionLocator;
} 

FQUALIFIER
G4bool GPPropagatorInField_IntersectChord( GPPropagatorInField *This,
					   GPThreeVector StartPointA, 
					   GPThreeVector EndPointB,
					   G4double      &NewSafety,
					   G4double      &LinearStepLength,
					   GPThreeVector &IntersectionPoint)
{
  // Calculate the direction and length of the chord AB
  //
  //********************************************************************
  //@@@G4FWP GPVIntersectionLocator is the base class of GPMultiLevelLocator
  //  return GPMultiLevelLocator_IntersectChord(This->fIntersectionLocator,
  return GPVIntersectionLocator_IntersectChord(This->fIntersectionLocator,
  				StartPointA,EndPointB,NewSafety,
				This->fPreviousSafety,This->fPreviousSftOrigin,
				LinearStepLength,IntersectionPoint,NULL);
}





