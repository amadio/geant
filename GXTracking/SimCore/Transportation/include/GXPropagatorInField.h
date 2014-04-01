#ifndef GXPropagatorInField_hh 
#define GXPropagatorInField_hh  1

#include "GPTypeDef.h"
#include "GPThreeVector.h"

#include "GPNavigator.h"

#include "GXFieldTrack.h"
#include "GXChordFinder.h" 
#include "GXFieldManager.h"
#include "GXMultiLevelLocator.h"

struct GXPropagatorInField
{
  // ----------------------------------------------------------------------
  //  DATA Members
  // ----------------------------------------------------------------------
  
  //  ==================================================================
  //  INVARIANTS - Must not change during tracking
  
  //  ** PARAMETERS -----------
  G4int    fMax_loop_count;
  // Limit for the number of sub-steps taken in one call to ComputeStep
  G4bool   fUseSafetyForOptimisation;
  
  //  Thresholds for identifying "abnormal" cases - which cause looping
  G4int    fActionThreshold_NoZeroSteps;       //  Threshold # - above it act
  G4int    fSevereActionThreshold_NoZeroSteps; //  Threshold # to act harshly
  G4int    fAbandonThreshold_NoZeroSteps;      //  Threshold # to abandon
  G4double fZeroStepThreshold; 
  // Threshold *length* for counting of tiny or 'zero' steps 
  
  G4double fLargestAcceptableStep;
  // Maximum size of a step - for optimization (and to avoid problems)
  //  ** End of PARAMETERS -----
  
  G4double kCarTolerance;
  // Geometrical tolerance defining surface thickness
  
  G4bool   fAllocatedLocator;                    //  Book-keeping
  
  //  --------------------------------------------------------
  //  ** Dependent Objects - to which work is delegated 
  
  GXFieldManager         *fDetectorFieldMgr; 
  // The  Field Manager of the whole Detector.  (default)
  
  //  G4VIntersectionLocator *fIntersectionLocator;
  GXMultiLevelLocator *fIntersectionLocator;
  // Refines candidate intersection
  
  //  G4VCurvedTrajectoryFilter* fpTrajectoryFilter;
  // The filter encapsulates the algorithm which selects which
  // intermediate points should be stored in a trajectory. 
  // When it is NULL, no intermediate points will be stored.
  // Else PIF::ComputeStep must submit (all) intermediate
  // points it calculates, to this filter.  (jacek 04/11/2002)
  
  GPNavigator            *fNavigator;
  // Set externally - only by tracking / run manager
  //
  //  ** End of Dependent Objects ----------------------------
  
  //  End of INVARIANTS 
  //  ==================================================================
  
  //  STATE information
  //  -----------------
  GXFieldManager *fCurrentFieldMgr;
      // The  Field Manager of the current volume (may be the global)
  G4bool         fSetFieldMgr;  // Has it been set for the current step
  
  // Parameters of current step
  G4double       fCharge, fInitialMomentumModulus, fMass;
  G4double       fEpsilonStep;        // Relative accuracy of current Step
  GXFieldTrack   End_PointAndTangent; // End point storage
  G4bool         fParticleIsLooping;
  G4int          fNoZeroStep;         //  Count of zero Steps
  
  // State used for Optimisation
  G4double       fFull_CurveLen_of_LastAttempt; 
  G4double       fLast_ProposedStepLength; 
  // Previous step information -- for use in adjust step size
  GPThreeVector  fPreviousSftOrigin;
  G4double       fPreviousSafety; 
  // Last safety origin & value: for optimisation
  
  G4int          fVerboseLevel;
  // For debuging purposes
  
};

extern "C" {

FQUALIFIER
void GXPropagatorInField_Constructor( GXPropagatorInField *This, 
				      GPNavigator    *theNavigator, 
				      GXFieldManager *detectorFieldMgr,
				      GXMultiLevelLocator *vLocator );

FQUALIFIER
void GXPropagatorInField_RefreshIntersectionLocator(GXPropagatorInField *This);

FQUALIFIER
G4double
GXPropagatorInField_ComputeStep(GXPropagatorInField *This,
				GXFieldTrack&      pFieldTrack,
				G4double           CurrentProposedStepLength,
				G4double&          currentSafety,   // IN/OUT
				GPVPhysicalVolume* pPhysVol ) ;

FQUALIFIER
void GXPropagatorInField_ClearPropagatorState(GXPropagatorInField *This);

FQUALIFIER
GXFieldManager* 
GXPropagatorInField_FindAndSetFieldManager(GXPropagatorInField *This);

FQUALIFIER
GXChordFinder* GXPropagatorInField_GetChordFinder( GXPropagatorInField *This );

FQUALIFIER
void GXPropagatorInField_SetChargeMomentumMass( GXPropagatorInField *This,  
						G4double Charge,    // in e+ units
						G4double Momentum,  // in GeV/c 
						G4double Mass);      // in ? units

FQUALIFIER
GPThreeVector  GXPropagatorInField_EndPosition( GXPropagatorInField *This );

FQUALIFIER
G4double GXPropagatorInField_GetEpsilonStep( GXPropagatorInField *This );

FQUALIFIER
void     GXPropagatorInField_SetEpsilonStep( GXPropagatorInField *This,
					     G4double newEps );

FQUALIFIER
G4bool   GXPropagatorInField_IsParticleLooping( GXPropagatorInField *This );

FQUALIFIER
G4int    GXPropagatorInField_GetMaxLoopCount( GXPropagatorInField *This );

FQUALIFIER
void GXPropagatorInField_SetMaxLoopCount( GXPropagatorInField *This,
					  G4int new_max ) ;

FQUALIFIER
G4double GXPropagatorInField_GetDeltaIntersection( GXPropagatorInField *This );

FQUALIFIER
G4double GXPropagatorInField_GetDeltaOneStep( GXPropagatorInField *This );

FQUALIFIER
GXFieldTrack GXPropagatorInField_GetEndState( GXPropagatorInField *This );

FQUALIFIER 
G4double GXPropagatorInField_GetMinimumEpsilonStep( GXPropagatorInField *This );

FQUALIFIER 
void      GXPropagatorInField_SetMinimumEpsilonStep(GXPropagatorInField *This,
						    G4double newEpsMin );

FQUALIFIER 
G4double  GXPropagatorInField_GetMaximumEpsilonStep( GXPropagatorInField *This );

FQUALIFIER 
void      GXPropagatorInField_SetMaximumEpsilonStep(GXPropagatorInField *This,
						    G4double newEpsMax );

FQUALIFIER
void GXPropagatorInField_SetLargestAcceptableStep(GXPropagatorInField *This,
						  G4double newBigDist );

FQUALIFIER
G4double  GXPropagatorInField_GetLargestAcceptableStep( GXPropagatorInField *This );

FQUALIFIER
GXFieldManager*  GXPropagatorInField_GetCurrentFieldManager( GXPropagatorInField *This );

FQUALIFIER
void GXPropagatorInField_SetThresholdNoZeroStep( GXPropagatorInField *This,
						 G4int noAct,
						 G4int noHarsh,
						 G4int noAbandon );

FQUALIFIER
G4int GXPropagatorInField_GetThresholdNoZeroSteps( GXPropagatorInField *This,
						   G4int i );

FQUALIFIER 
G4double  GXPropagatorInField_GetZeroStepThreshold( GXPropagatorInField *This );

FQUALIFIER 
void      GXPropagatorInField_SetZeroStepThreshold(GXPropagatorInField *This,
						   G4double newLength );

FQUALIFIER
void GXPropagatorInField_SetDetectorFieldManager(GXPropagatorInField *This,
						 GXFieldManager* newDetectorFieldManager);

FQUALIFIER
void GXPropagatorInField_SetUseSafetyForOptimization(GXPropagatorInField *This,
						     G4bool value );

FQUALIFIER 
G4bool 
GXPropagatorInField_GetUseSafetyForOptimization( GXPropagatorInField *This ); 

FQUALIFIER 
void GXPropagatorInField_SetNavigatorForPropagating(GXPropagatorInField *This, 
						    GPNavigator *SimpleOrMultiNavigator );

FQUALIFIER
GPNavigator* GXPropagatorInField_GetNavigatorForPropagating( GXPropagatorInField *This );

FQUALIFIER 
void GXPropagatorInField_SetIntersectionLocator(GXPropagatorInField *This,
						GXMultiLevelLocator *pIntLoc );

FQUALIFIER
GXMultiLevelLocator* GXPropagatorInField_GetIntersectionLocator( GXPropagatorInField *This );

FQUALIFIER
G4bool GXPropagatorInField_IntersectChord( GXPropagatorInField *This,
					   GPThreeVector StartPointA, 
					   GPThreeVector EndPointB,
					   G4double      &NewSafety,
					   G4double      &LinearStepLength,
					   GPThreeVector &IntersectionPoint);

}

#endif 
