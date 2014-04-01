#ifndef GXMULTILEVELLOCATOR_HH
#define GXMULTILEVELLOCATOR_HH

//dummy struct for now

#include "GPGeomdefs.h" 
#include "GXChordFinder.h"
#include "GPNavigator.h"
#include "GPTouchableHistory.h"

struct GXMultiLevelLocator
{
  //from base G4VIntersectionLocator

  //G4double kCarTolerance;         // Constant
  //G4int    fVerboseLevel;          // For debugging
  G4bool   fUseNormalCorrection;   // Configuration parameter

  G4double       fiEpsilonStep;
  G4double       fiDeltaIntersection;
  G4bool         fiUseSafety;

  GXChordFinder *fiChordFinder;
  GPNavigator   *fiNavigator;

  //Parameters at each physical step by calling method  by G4PropagatorInField
  GPNavigator *fHelpingNavigator; // Helper for location
  GPTouchableHistory *fpTouchable; // Touchable history hook

  //addition for G4MultiLevelLocator
  //  static const G4int max_depth=10; //move to GPGeomdefs.h
  GXFieldTrack* ptrInterMedFT[max_depth+1];
};

extern "C" {

FQUALIFIER
void GXMultiLevelLocator_Constructor( GXMultiLevelLocator *This,
				      GPNavigator *theNavigator);

FQUALIFIER
void GPVIntersectionLocator_GPVIntersectionLocator(GXMultiLevelLocator *This,
					      GPNavigator *theNavigator);

FQUALIFIER
G4bool GXMultiLevelLocator_EstimateIntersectionPoint(GXMultiLevelLocator *This, 
		GXFieldTrack&       CurveStartPointVelocity,       // A
                GXFieldTrack&       CurveEndPointVelocity,         // B
                GPThreeVector&      TrialPoint,                    // E
                GXFieldTrack&       IntersectedOrRecalculatedFT,   // Output
                G4bool&             recalculatedEndPoint,          // Out
                G4double&           previousSafety,                // In/Out
		GPThreeVector&      previousSftOrigin);             // In/Out

FQUALIFIER
GXFieldTrack 
GPVIntersectionLocator_ReEstimateEndpoint(GXMultiLevelLocator *This, 
					  GXFieldTrack &CurrentStateA,  
					  GXFieldTrack &EstimatedEndStateB,
					  G4double      linearDistSq,
					  G4double      curveDist );
FQUALIFIER
GPThreeVector GPVIntersectionLocator_GetLocalSurfaceNormal(GXMultiLevelLocator *This, 
		      GPThreeVector &CurrentE_Point, G4bool &validNormal);

FQUALIFIER
G4bool GPVIntersectionLocator_AdjustmentOfFoundIntersection( GXMultiLevelLocator *This, 
			       GPThreeVector &CurrentA_Point,
                               GPThreeVector &CurrentE_Point, 
                               GPThreeVector &CurrentF_Point,
                               GPThreeVector &MomentumDir,
                               G4bool         IntersectAF,
			       GPThreeVector &IntersectionPoint,  // I/O
			       G4double      &NewSafety,          // I/O 
			       G4double      &fPreviousSafety,    // I/O
			       GPThreeVector &fPreviousSftOrigin );// I/O
FQUALIFIER
GPThreeVector GPVIntersectionLocator_GetSurfaceNormal(GXMultiLevelLocator *This,
				     GPThreeVector &CurrentInt_Point,
				     G4bool &validNormal);

FQUALIFIER
GPThreeVector GPVIntersectionLocator_GetGlobalSurfaceNormal(
				     GXMultiLevelLocator *This, 
		                     GPThreeVector &CurrentE_Point,
		                     G4bool &validNormal);

FQUALIFIER
GPThreeVector 
GPVIntersectionLocator_GetLastSurfaceNormal( GXMultiLevelLocator *This, 
					     GPThreeVector intersectPoint,
					     G4bool &normalIsValid);

FQUALIFIER 
G4double GPVIntersectionLocator_GetDeltaIntersectionFor(GXMultiLevelLocator *This);

FQUALIFIER 
G4double GPVIntersectionLocator_GetEpsilonStepFor(GXMultiLevelLocator *This);

FQUALIFIER 
GPNavigator* GPVIntersectionLocator_GetNavigatorFor(GXMultiLevelLocator *This);

FQUALIFIER
GXChordFinder* GPVIntersectionLocator_GetChordFinderFor(GXMultiLevelLocator *This);

FQUALIFIER 
G4int GPVIntersectionLocator_GetVerboseFor(GXMultiLevelLocator *This);

FQUALIFIER 
G4bool GPVIntersectionLocator_GetAdjustementOfFoundIntersection(GXMultiLevelLocator *This );

FQUALIFIER 
void GPVIntersectionLocator_AddAdjustementOfFoundIntersection(
				  GXMultiLevelLocator *This, 
				  G4bool UseCorrection );

FQUALIFIER 
void GPVIntersectionLocator_SetEpsilonStepFor(GXMultiLevelLocator *This,
					      G4double EpsilonStep );

FQUALIFIER 
void GPVIntersectionLocator_SetDeltaIntersectionFor(GXMultiLevelLocator *This, 
						    G4double deltaIntersection );

FQUALIFIER 
void GPVIntersectionLocator_SetNavigatorFor( GXMultiLevelLocator *This, 
					     GPNavigator *fNavigator );

FQUALIFIER 
void GPVIntersectionLocator_SetChordFinderFor(GXMultiLevelLocator *This, 
					      GXChordFinder *fCFinder );

FQUALIFIER 
void GPVIntersectionLocator_SetSafetyParametersFor(GXMultiLevelLocator *This,
						   G4bool UseSafety );

FQUALIFIER 
void GPVIntersectionLocator_SetVerboseFor(GXMultiLevelLocator *This, 
					  G4int fVerbose);

FQUALIFIER G4bool
GPVIntersectionLocator_IntersectChord( GXMultiLevelLocator *This, 
				       GPThreeVector  StartPointA, 
				       GPThreeVector  EndPointB,
				       G4double      &NewSafety,
				       G4double      &PreviousSafety,
				       GPThreeVector &PreviousSftOrigin,
				       G4double      &LinearStepLength,
				       GPThreeVector &IntersectionPoint,
				       G4bool        *ptrCalledNavigator );


}

#endif
