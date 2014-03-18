#include "GXPropagatorInField.h"
#include "GXChordFinder.h"
#include "GPConstants.h"
#include "GPUtils.h"

///////////////////////////////////////////////////////////////////////////
//
// Constructors and destructor

FQUALIFIER
void GXPropagatorInField_Constructor( GXPropagatorInField *This, 
				      GPNavigator    *theNavigator, 
				      GXFieldManager *detectorFieldMgr,
				      GXMultiLevelLocator *vLocator )
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

  GXFieldTrack_Constructor(&(This->End_PointAndTangent),
			   GPThreeVector_create(0.,0.,0.),
			   GPThreeVector_create(0.,0.,0.),
			   0.0,0.0,0.0,0.0);

  This->fParticleIsLooping = false;
  This->fNoZeroStep = 0;
  This->fVerboseLevel = 0;

  if(This->fDetectorFieldMgr) 
  { 
    This->fEpsilonStep = 
      GXFieldManager_GetMaximumEpsilonStep(This->fDetectorFieldMgr);
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
    GXMultiLevelLocator mLocator;
    GXMultiLevelLocator_Constructor(&mLocator,theNavigator);
    This->fIntersectionLocator = &mLocator;
    This->fAllocatedLocator=true;
  }else{
    This->fIntersectionLocator=vLocator;
    This->fAllocatedLocator=false;
  }

  GXPropagatorInField_RefreshIntersectionLocator(This);//Copy all relevant parameters 
}

//********************************************************************
// Update the IntersectionLocator with current parameters
FQUALIFIER
void GXPropagatorInField_RefreshIntersectionLocator(GXPropagatorInField *This)
{
  GPVIntersectionLocator_SetEpsilonStepFor(This->fIntersectionLocator,
					   This->fEpsilonStep);
  GPVIntersectionLocator_SetDeltaIntersectionFor(This->fIntersectionLocator,
		 GXFieldManager_GetDeltaIntersection(This->fCurrentFieldMgr));
  GPVIntersectionLocator_SetChordFinderFor(This->fIntersectionLocator,
			 GXPropagatorInField_GetChordFinder(This));
  GPVIntersectionLocator_SetSafetyParametersFor(This->fIntersectionLocator, 
					       This->fUseSafetyForOptimisation);
}
///////////////////////////////////////////////////////////////////////////
//
// Compute the next geometric Step

FQUALIFIER
G4double
GXPropagatorInField_ComputeStep(GXPropagatorInField *This,
				GXFieldTrack&      pFieldTrack,
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
    This->fCurrentFieldMgr = GXPropagatorInField_FindAndSetFieldManager(This); 
  //@@@G4FWP*****************************************************************
  // For the next call, the field manager must again be set
  This->fSetFieldMgr= false;

  // GetChordFinder()->SetChargeMomentumMass(fCharge, 
  //                                         fInitialMomentumModulus, fMass); 
  GXChordFinder *chordFinder = 
    GXFieldManager_GetChordFinder(This->fCurrentFieldMgr);

  /*
  GXChordFinder_SetChargeMomentumMass(chordFinder,
				      This->fCharge, 
				      This->fInitialMomentumModulus, 
  				      This->fMass);
  */

  // Values for Intersection Locator has to be updated on each call for the
  // case that CurrentFieldManager has changed from the one of previous step
  //  GXPropagatorInField_RefreshIntersectionLocator(This);

  //  G4FieldTrack  CurrentState(pFieldTrack);
  GXFieldTrack&  CurrentState  = pFieldTrack;
  //  G4FieldTrack  OriginalState = CurrentState;
  GXFieldTrack&  OriginalState = pFieldTrack;

  // If the Step length is "infinite", then an approximate-maximum Step
  // length (used to calculate the relative accuracy) must be guessed.
  //
  if( CurrentProposedStepLength >= This->fLargestAcceptableStep )
  {
    GPThreeVector StartPointA, VelocityUnit;
    StartPointA  = GXFieldTrack_GetPosition(&pFieldTrack);
    //    VelocityUnit = GXFieldTrack_GetMomentumDir(&pFieldTrack);
    VelocityUnit = GPThreeVector_unit(GXFieldTrack_GetMomentum(&pFieldTrack));

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
    GXFieldManager_GetDeltaOneStep(This->fCurrentFieldMgr)/CurrentProposedStepLength;

  // G4double raw_epsilon= epsilon;
  G4double epsilonMin= GXFieldManager_GetMinimumEpsilonStep(This->fCurrentFieldMgr);
  G4double epsilonMax= GXFieldManager_GetMaximumEpsilonStep(This->fCurrentFieldMgr);
  if( epsilon < epsilonMin ) epsilon = epsilonMin;
  if( epsilon > epsilonMax ) epsilon = epsilonMax;
  GXPropagatorInField_SetEpsilonStep(This, epsilon );

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
       //G4Exception("GXPropagatorInField_ComputeStep()", "GeomNav1002",
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
    GXFieldTrack& SubStepStartState = CurrentState;

    GPThreeVector SubStartPoint = GXFieldTrack_GetPosition(&CurrentState); 

    if( !first_substep) {
      //      GPNavigator_LocateGlobalPointWithinVolume(This->fNavigator,SubStartPoint);
    }

    // How far to attempt to move the particle !
    //
    h_TrialStepSize = CurrentProposedStepLength - StepTaken;

    // Integrate as far as "chord miss" rule allows.
    //

    s_length_taken = GXChordFinder_AdvanceChordLimited(
                             GXPropagatorInField_GetChordFinder(This), 
                             CurrentState,    // Position & velocity
                             h_TrialStepSize,
                             This->fEpsilonStep);
    //                             This->fPreviousSftOrigin,
    //                             This->fPreviousSafety
    //                             );
    //  CurrentState is now updated with the final position and velocity. 

    This->fFull_CurveLen_of_LastAttempt = s_length_taken;

    GPThreeVector  EndPointB = GXFieldTrack_GetPosition(&CurrentState); 
    GPThreeVector  InterSectionPointE;
    G4double       LinearStepLength;
 
    // Intersect chord AB with geometry
    intersects= GXPropagatorInField_IntersectChord(This, 
    						   SubStartPoint, EndPointB, 
    						   NewSafety, LinearStepLength, 
    						   InterSectionPointE );
    // E <- Intersection Point of chord AB and either volume A's surface 
    //                                  or a daughter volume's surface ..

    if( first_substep ) { 
       currentSafety = NewSafety;
    } // Updating safety in other steps is potential future extention

    if( intersects )
    {
      GXFieldTrack& IntersectPointVelct_G = CurrentState;  // FT-Def-Construct

      // Find the intersection point of AB true path with the surface
      //   of vol(A), if it exists. Start with point E as first "estimate".
      G4bool recalculatedEndPt= false;

      G4bool found_intersection = false;
	GXMultiLevelLocator_EstimateIntersectionPoint(
                                    This->fIntersectionLocator,
				    SubStepStartState, CurrentState, 
				    InterSectionPointE, IntersectPointVelct_G,
				    recalculatedEndPt,This->fPreviousSafety,
				    This->fPreviousSftOrigin);
      intersects = intersects && found_intersection;
      if( found_intersection ) {        
	This->End_PointAndTangent= IntersectPointVelct_G;  //G is our EndPoint
	StepTaken = TruePathLength = 
	    GXFieldTrack_GetCurveLength(&IntersectPointVelct_G)
	  - GXFieldTrack_GetCurveLength(&OriginalState);
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
     //G4Exception("GXPropagatorInField_ComputeStep()",
     //            "GeomNav1002", JustWarning, message);
     This->fNoZeroStep = 0; 
  }
 
  return TruePathLength;
}

FQUALIFIER
void GXPropagatorInField_ClearPropagatorState(GXPropagatorInField *This)
{
  // Goal: Clear all memory of previous steps,  cached information

  This->fParticleIsLooping= false;
  This->fNoZeroStep= 0;

  GXFieldTrack_Constructor(&(This->End_PointAndTangent), 
			    GPThreeVector_create(0.,0.,0.),
			    GPThreeVector_create(0.,0.,0.),
			    0.0,0.0,0.0,0.0); 
  This->fFull_CurveLen_of_LastAttempt = -1; 
  This->fLast_ProposedStepLength = -1;

  This->fPreviousSftOrigin= GPThreeVector_create(0.,0.,0.);
  This->fPreviousSafety= 0.0;
}

//G4FieldManager* 
//GXPropagatorInField_FindAndSetFieldManager( G4VPhysicalVolume* pCurrentPhysicalVolume)

FQUALIFIER
GXFieldManager* 
GXPropagatorInField_FindAndSetFieldManager(GXPropagatorInField *This)
{
  GXFieldManager* currentFieldMgr;

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

// ------------------------------------------------------------------------
// G4PropagatorInField.icc inline implementation

FQUALIFIER
GXChordFinder* GXPropagatorInField_GetChordFinder( GXPropagatorInField *This )
{
  // The "Chord Finder" of the current Field Mgr is used
  //    -- this could be of the global field manager
  //        or that of another, from the current volume 
  return GXFieldManager_GetChordFinder( This->fCurrentFieldMgr ); 
}

FQUALIFIER
void GXPropagatorInField_SetChargeMomentumMass( GXPropagatorInField *This,  
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
GPThreeVector  GXPropagatorInField_EndPosition( GXPropagatorInField *This )
{
  return  GXFieldTrack_GetPosition(&(This->End_PointAndTangent)); 
}

FQUALIFIER
G4double GXPropagatorInField_GetEpsilonStep( GXPropagatorInField *This )
{ 
  return This->fEpsilonStep; 
}

FQUALIFIER
void     GXPropagatorInField_SetEpsilonStep( GXPropagatorInField *This,
					     G4double newEps )
{
  This->fEpsilonStep=newEps;
}

FQUALIFIER
G4bool   GXPropagatorInField_IsParticleLooping( GXPropagatorInField *This )
{
  return This->fParticleIsLooping;
}

FQUALIFIER
G4int    GXPropagatorInField_GetMaxLoopCount( GXPropagatorInField *This )
{
  return This->fMax_loop_count;
}

FQUALIFIER
void GXPropagatorInField_SetMaxLoopCount( GXPropagatorInField *This,
					      G4int new_max ) 
{
  This->fMax_loop_count = new_max;
}

// #if 0
FQUALIFIER
G4double GXPropagatorInField_GetDeltaIntersection( GXPropagatorInField *This )
{
  return GXFieldManager_GetDeltaIntersection(This->fCurrentFieldMgr);
} 

FQUALIFIER
G4double GXPropagatorInField_GetDeltaOneStep( GXPropagatorInField *This )
{
  return GXFieldManager_GetDeltaOneStep(This->fCurrentFieldMgr);
}

FQUALIFIER
GXFieldTrack GXPropagatorInField_GetEndState( GXPropagatorInField *This )
{
  return This->End_PointAndTangent;
}

// Minimum for Relative accuracy of a Step in volumes of global field
FQUALIFIER 
G4double GXPropagatorInField_GetMinimumEpsilonStep( GXPropagatorInField *This )
{
  return GXFieldManager_GetMinimumEpsilonStep(This->fDetectorFieldMgr);
}

FQUALIFIER 
void      GXPropagatorInField_SetMinimumEpsilonStep(GXPropagatorInField *This,
						    G4double newEpsMin )
{
  GXFieldManager_SetMinimumEpsilonStep(This->fDetectorFieldMgr, newEpsMin);
}

// Maximum for Relative accuracy of any Step 
FQUALIFIER 
G4double  GXPropagatorInField_GetMaximumEpsilonStep( GXPropagatorInField *This )
{
  return GXFieldManager_GetMaximumEpsilonStep(This->fDetectorFieldMgr);
}

FQUALIFIER 
void      GXPropagatorInField_SetMaximumEpsilonStep(GXPropagatorInField *This,
						    G4double newEpsMax )
{
  GXFieldManager_SetMaximumEpsilonStep(This->fDetectorFieldMgr, newEpsMax );
}

FQUALIFIER
void GXPropagatorInField_SetLargestAcceptableStep(GXPropagatorInField *This,
						  G4double newBigDist )
{
  if( This->fLargestAcceptableStep>0.0 )
  {
    This->fLargestAcceptableStep = newBigDist;
  }
}

FQUALIFIER
G4double  GXPropagatorInField_GetLargestAcceptableStep( GXPropagatorInField *This )
{
  return This->fLargestAcceptableStep;
}

FQUALIFIER
GXFieldManager*  GXPropagatorInField_GetCurrentFieldManager( GXPropagatorInField *This )
{
  return This->fCurrentFieldMgr;
} 
FQUALIFIER
void GXPropagatorInField_SetThresholdNoZeroStep( GXPropagatorInField *This,
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
G4int GXPropagatorInField_GetThresholdNoZeroSteps( GXPropagatorInField *This,
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
G4double  GXPropagatorInField_GetZeroStepThreshold( GXPropagatorInField *This )
{ 
  return This->fZeroStepThreshold; 
}

FQUALIFIER 
void      GXPropagatorInField_SetZeroStepThreshold(GXPropagatorInField *This,
						   G4double newLength )
{ 
  This->fZeroStepThreshold= newLength;
}

FQUALIFIER
void GXPropagatorInField_SetDetectorFieldManager(GXPropagatorInField *This,
				      GXFieldManager* newDetectorFieldManager)
{
  This->fDetectorFieldMgr = newDetectorFieldManager; 
}

FQUALIFIER
void GXPropagatorInField_SetUseSafetyForOptimization(GXPropagatorInField *This,
						     G4bool value )
{
  This->fUseSafetyForOptimisation= value;
}

FQUALIFIER 
G4bool 
GXPropagatorInField_GetUseSafetyForOptimization( GXPropagatorInField *This ) 
{ 
  return This->fUseSafetyForOptimisation; 
}

FQUALIFIER 
void GXPropagatorInField_SetNavigatorForPropagating(GXPropagatorInField *This, 
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
GPNavigator* GXPropagatorInField_GetNavigatorForPropagating( GXPropagatorInField *This )
{
  return This->fNavigator;
} 

FQUALIFIER 
void GXPropagatorInField_SetIntersectionLocator(GXPropagatorInField *This,
						GXMultiLevelLocator *pIntLoc )
{
  if(pIntLoc)  { This->fIntersectionLocator= pIntLoc; }
}

//G4VIntersectionLocator* GXPropagatorInField_GetIntersectionLocator( GXPropagatorInField *This )
FQUALIFIER
GXMultiLevelLocator* GXPropagatorInField_GetIntersectionLocator( GXPropagatorInField *This )
{
  return This->fIntersectionLocator;
} 

FQUALIFIER
G4bool GXPropagatorInField_IntersectChord( GXPropagatorInField *This,
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
