#ifndef GXCHORDFINDER_HH
#define GXCHORDFINDER_HH

#include "GXMagIntegratorDriver.h"
#include "GXClassicalRK4.h"
#include "GXEquationOfMotion.h"
#include "GXFieldTrack.h"
#include "GXChordFinder.h"

struct GXChordFinder
{ 
  G4double fDeltaChord; //  Maximum miss distance 

  G4double fFirstFraction; // Internal parameters
  G4double fFractionLast;
  G4double fFractionNextEstimate;
  G4double fLastStepEstimate_Unconstrained;
  
  //  DEPENDENT Objects
  GXEquationOfMotion* fEquation; 
  GXClassicalRK4*     fDriversStepper; 
  GXMagInt_Driver*    fIntgrDriver;
};

extern "C" {

FQUALIFIER
void GXChordFinder_Constructor(GXChordFinder *This,
			       GXMagInt_Driver* pIntegrationDriver);

FQUALIFIER
void GXChordFinder_Constructor2( GXChordFinder *This,
				 GXMagneticField*        theMagField,
				 G4double                stepMinimum, 
				 GXClassicalRK4* pItsStepper );

FQUALIFIER
G4double 
GXChordFinder_AdvanceChordLimited2( GXChordFinder *This,
				   GXFieldTrack& yCurrent,
				   G4double      stepMax,
				   G4double      epsStep);

FQUALIFIER
G4double 
GXChordFinder_AdvanceChordLimited( GXChordFinder *This,
				   GXFieldTrack& yCurrent,
				   G4double      stepMax,
				   G4double      epsStep);

FQUALIFIER
G4double
GXChordFinder_FindNextChord( GXChordFinder *This, 
			     const  GXFieldTrack& yStart,
			     G4double     stepMax,
			     GXFieldTrack&   yEnd, // Endpoint
			     G4double&   dyErrPos, // Error of endpoint
			     G4double    epsStep,
			     G4double*  pStepForAccuracy);

FQUALIFIER
G4double
GXChordFinder_FindNextChord2( GXChordFinder *This, 
			     const  GXFieldTrack& yStart,
			     G4double     stepMax,
			     GXFieldTrack&   yEnd, // Endpoint
			     G4double&   dyErrPos, // Error of endpoint
			     G4double    epsStep,
			     G4double*  pStepForAccuracy);


FQUALIFIER
G4double GXChordFinder_NewStep( GXChordFinder *This,
				G4double  stepTrialOld, 
                                G4double  dChordStep, // Curr. dchord achieved
                                G4double& stepEstimate_Unconstrained ) ;

FQUALIFIER
GXFieldTrack
GXChordFinder_ApproxCurvePointS( GXChordFinder *This, 
				 const GXFieldTrack&  CurveA_PointVelocity, 
				 const GXFieldTrack&  CurveB_PointVelocity, 
				 const GXFieldTrack&  ApproxCurveV,
				 const GPThreeVector& CurrentE_Point,
				 const GPThreeVector& CurrentF_Point,
				 const GPThreeVector& PointG,
				 G4bool first, G4double eps_step);

FQUALIFIER
GXFieldTrack 
GXChordFinder_ApproxCurvePointV( GXChordFinder *This, 
				 const GXFieldTrack& CurveA_PointVelocity, 
				 const GXFieldTrack& CurveB_PointVelocity, 
				 const GPThreeVector& CurrentE_Point,
				 G4double eps_step);

FQUALIFIER
GXMagInt_Driver* GXChordFinder_GetIntegrationDriver(GXChordFinder *This);

FQUALIFIER
G4bool GXChordFinder_AcceptableMissDist(GXChordFinder *This,
					G4double dChordStep);

FQUALIFIER
G4double GXChordFinder_GetDeltaChord(GXChordFinder *This);

FQUALIFIER
void GXChordFinder_SetDeltaChord(GXChordFinder *This,G4double newval);

FQUALIFIER 
G4double GXChordFinder_GetLastStepEstimateUnc(GXChordFinder *This);

FQUALIFIER 
void GXChordFinder_SetLastStepEstimateUnc(GXChordFinder *This,
					  G4double stepEst );

FQUALIFIER
void GXChordFinder_ResetStepEstimate(GXChordFinder *This);


FQUALIFIER 
G4double GXChordFinder_GetFirstFraction(GXChordFinder *This);

FQUALIFIER 
G4double GXChordFinder_GetFractionLast(GXChordFinder *This);

FQUALIFIER 
G4double GXChordFinder_InvParabolic ( const G4double xa, 
				      const G4double ya,
				      const G4double xb, 
				      const G4double yb,
				      const G4double xc, 
				      const G4double yc );

}

#endif
