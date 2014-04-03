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
// $Id: G4MultiLevelLocator.hh,v 1.3 2010-07-13 15:59:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Class G4MultiLevelLocator 
//
// class description:
// 
// Implementing the calculation of the intersection point with a boundary when
// PropagationInField is used. Derived from method LocateIntersectionPoint()
// from G4PropagatorInField, it is based on a linear method for finding the
// intersection point by means of a 'depth' algorithm in case of slow progress
// (intersection is not found after 100 trials).

// History:
// -------
// 27.10.08 - Tatiana Nikitina: Derived from LocateIntersectionPoint() from 
//                              G4PropagatorInField class
// ---------------------------------------------------------------------------

#ifndef GPMULTILEVELLOCATOR_HH
#define GPMULTILEVELLOCATOR_HH

//dummy struct for now

#include "GPGeomdefs.h" 
#include "GPChordFinder.h"
#include "GPNavigator.h"
#include "GPTouchableHistory.h"

struct GPMultiLevelLocator
{
  //from base G4VIntersectionLocator

  //G4double kCarTolerance;         // Constant
  //G4int    fVerboseLevel;          // For debugging
  G4bool   fUseNormalCorrection;   // Configuration parameter

  G4double       fiEpsilonStep;
  G4double       fiDeltaIntersection;
  G4bool         fiUseSafety;

  GPChordFinder *fiChordFinder;
  GPNavigator   *fiNavigator;

  //Parameters at each physical step by calling method  by G4PropagatorInField
  GPNavigator *fHelpingNavigator; // Helper for location
  GPTouchableHistory *fpTouchable; // Touchable history hook

  //addition for G4MultiLevelLocator
  //  static const G4int max_depth=10; //move to GPGeomdefs.h
  GPFieldTrack* ptrInterMedFT[max_depth+1];
};

extern "C" {

FQUALIFIER
void GPMultiLevelLocator_Constructor( GPMultiLevelLocator *This,
				      GPNavigator *theNavigator);

FQUALIFIER
void GPVIntersectionLocator_GPVIntersectionLocator(GPMultiLevelLocator *This,
					      GPNavigator *theNavigator);

FQUALIFIER
G4bool GPMultiLevelLocator_EstimateIntersectionPoint(GPMultiLevelLocator *This, 
		GPFieldTrack&       CurveStartPointVelocity,       // A
                GPFieldTrack&       CurveEndPointVelocity,         // B
                GPThreeVector&      TrialPoint,                    // E
                GPFieldTrack&       IntersectedOrRecalculatedFT,   // Output
                G4bool&             recalculatedEndPoint,          // Out
                G4double&           previousSafety,                // In/Out
		GPThreeVector&      previousSftOrigin);             // In/Out

FQUALIFIER
GPFieldTrack 
GPVIntersectionLocator_ReEstimateEndpoint(GPMultiLevelLocator *This, 
					  GPFieldTrack &CurrentStateA,  
					  GPFieldTrack &EstimatedEndStateB,
					  G4double      linearDistSq,
					  G4double      curveDist );
FQUALIFIER
GPThreeVector GPVIntersectionLocator_GetLocalSurfaceNormal(GPMultiLevelLocator *This, 
		      GPThreeVector &CurrentE_Point, G4bool &validNormal);

FQUALIFIER
G4bool GPVIntersectionLocator_AdjustmentOfFoundIntersection( GPMultiLevelLocator *This, 
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
GPThreeVector GPVIntersectionLocator_GetSurfaceNormal(GPMultiLevelLocator *This,
				     GPThreeVector &CurrentInt_Point,
				     G4bool &validNormal);

FQUALIFIER
GPThreeVector GPVIntersectionLocator_GetGlobalSurfaceNormal(
				     GPMultiLevelLocator *This, 
		                     GPThreeVector &CurrentE_Point,
		                     G4bool &validNormal);

FQUALIFIER
GPThreeVector 
GPVIntersectionLocator_GetLastSurfaceNormal( GPMultiLevelLocator *This, 
					     GPThreeVector intersectPoint,
					     G4bool &normalIsValid);

FQUALIFIER 
G4double GPVIntersectionLocator_GetDeltaIntersectionFor(GPMultiLevelLocator *This);

FQUALIFIER 
G4double GPVIntersectionLocator_GetEpsilonStepFor(GPMultiLevelLocator *This);

FQUALIFIER 
GPNavigator* GPVIntersectionLocator_GetNavigatorFor(GPMultiLevelLocator *This);

FQUALIFIER
GPChordFinder* GPVIntersectionLocator_GetChordFinderFor(GPMultiLevelLocator *This);

FQUALIFIER 
G4int GPVIntersectionLocator_GetVerboseFor(GPMultiLevelLocator *This);

FQUALIFIER 
G4bool GPVIntersectionLocator_GetAdjustementOfFoundIntersection(GPMultiLevelLocator *This );

FQUALIFIER 
void GPVIntersectionLocator_AddAdjustementOfFoundIntersection(
				  GPMultiLevelLocator *This, 
				  G4bool UseCorrection );

FQUALIFIER 
void GPVIntersectionLocator_SetEpsilonStepFor(GPMultiLevelLocator *This,
					      G4double EpsilonStep );

FQUALIFIER 
void GPVIntersectionLocator_SetDeltaIntersectionFor(GPMultiLevelLocator *This, 
						    G4double deltaIntersection );

FQUALIFIER 
void GPVIntersectionLocator_SetNavigatorFor( GPMultiLevelLocator *This, 
					     GPNavigator *fNavigator );

FQUALIFIER 
void GPVIntersectionLocator_SetChordFinderFor(GPMultiLevelLocator *This, 
					      GPChordFinder *fCFinder );

FQUALIFIER 
void GPVIntersectionLocator_SetSafetyParametersFor(GPMultiLevelLocator *This,
						   G4bool UseSafety );

FQUALIFIER 
void GPVIntersectionLocator_SetVerboseFor(GPMultiLevelLocator *This, 
					  G4int fVerbose);

FQUALIFIER G4bool
GPVIntersectionLocator_IntersectChord( GPMultiLevelLocator *This, 
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
