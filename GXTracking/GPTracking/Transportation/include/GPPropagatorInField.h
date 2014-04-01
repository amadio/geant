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
// $Id: G4PropagatorInField.hh,v 1.19 2009-11-13 17:34:26 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Class G4PropagatorInField 
//
// class description:
// 
// This class performs the navigation/propagation of a particle/track 
// in a magnetic field. The field is in general non-uniform.
// For the calculation of the path, it relies on the class G4ChordFinder.
//
// Key Method: ComputeStep(..)

// History:
// -------
// 25.10.96 John Apostolakis,  design and implementation 
// 25.03.97 John Apostolakis,  adaptation for G4Transportation and cleanup
//  8.11.02 John Apostolakis,  changes to enable use of safety in intersecting
// ---------------------------------------------------------------------------

#ifndef GPPropagatorInField_hh 
#define GPPropagatorInField_hh  1

//#include "G4Types.hh"
#include "GPTypeDef.h"
#include "GPThreeVector.h"

//#include <vector>

//#include "G4FieldTrack.hh"
//#include "G4FieldManager.hh"
//#include "G4VIntersectionLocator.hh"
#include "GPFieldTrack.h"
#include "GPFieldManager.h"
#include "GPMultiLevelLocator.h"
#include "GPNavigator.h"
#include "GPChordFinder.h" 

//class G4ChordFinder; 

//class G4Navigator;
//class G4VPhysicalVolume;
//class G4VCurvedTrajectoryFilter;

struct GPPropagatorInField
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
  
  GPFieldManager         *fDetectorFieldMgr; 
  // The  Field Manager of the whole Detector.  (default)
  
  //  G4VIntersectionLocator *fIntersectionLocator;
  GPMultiLevelLocator *fIntersectionLocator;
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
  GPFieldManager *fCurrentFieldMgr;
      // The  Field Manager of the current volume (may be the global)
  G4bool         fSetFieldMgr;  // Has it been set for the current step
  
  // Parameters of current step
  G4double       fCharge, fInitialMomentumModulus, fMass;
  G4double       fEpsilonStep;        // Relative accuracy of current Step
  GPFieldTrack   End_PointAndTangent; // End point storage
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
void GPPropagatorInField_Constructor( GPPropagatorInField *This, 
				      GPNavigator    *theNavigator, 
				      GPFieldManager *detectorFieldMgr,
				      GPMultiLevelLocator *vLocator );

FQUALIFIER
void GPPropagatorInField_RefreshIntersectionLocator(GPPropagatorInField *This);

FQUALIFIER
G4double
GPPropagatorInField_ComputeStep(GPPropagatorInField *This,
				GPFieldTrack&      pFieldTrack,
				G4double           CurrentProposedStepLength,
				G4double&          currentSafety,   // IN/OUT
				GPVPhysicalVolume* pPhysVol ) ;

FQUALIFIER
void GPPropagatorInField_ClearPropagatorState(GPPropagatorInField *This);

FQUALIFIER
GPFieldManager* 
GPPropagatorInField_FindAndSetFieldManager(GPPropagatorInField *This);

FQUALIFIER
G4int GPPropagatorInField_SetVerboseLevel(GPPropagatorInField *This, G4int level );

FQUALIFIER
GPChordFinder* GPPropagatorInField_GetChordFinder( GPPropagatorInField *This );

FQUALIFIER
void GPPropagatorInField_SetChargeMomentumMass( GPPropagatorInField *This,  
						G4double Charge,    // in e+ units
						G4double Momentum,  // in GeV/c 
						G4double Mass);      // in ? units

FQUALIFIER
GPThreeVector  GPPropagatorInField_EndPosition( GPPropagatorInField *This );

FQUALIFIER
GPThreeVector  GPPropagatorInField_EndMomentumDir( GPPropagatorInField *This );

FQUALIFIER
G4double GPPropagatorInField_GetEpsilonStep( GPPropagatorInField *This );

FQUALIFIER
void     GPPropagatorInField_SetEpsilonStep( GPPropagatorInField *This,
					     G4double newEps );

FQUALIFIER
G4bool   GPPropagatorInField_IsParticleLooping( GPPropagatorInField *This );

FQUALIFIER
G4int    GPPropagatorInField_GetMaxLoopCount( GPPropagatorInField *This );

FQUALIFIER
void GPPropagatorInField_SetMaxLoopCount( GPPropagatorInField *This,
					  G4int new_max ) ;

FQUALIFIER
G4double GPPropagatorInField_GetDeltaIntersection( GPPropagatorInField *This );

FQUALIFIER
G4double GPPropagatorInField_GetDeltaOneStep( GPPropagatorInField *This );

FQUALIFIER
G4int GPPropagatorInField_GetVerboseLevel( GPPropagatorInField *This );

FQUALIFIER
G4int GPPropagatorInField_Verbose( GPPropagatorInField *This );

FQUALIFIER
GPFieldTrack GPPropagatorInField_GetEndState( GPPropagatorInField *This );

FQUALIFIER 
G4double GPPropagatorInField_GetMinimumEpsilonStep( GPPropagatorInField *This );

FQUALIFIER 
void      GPPropagatorInField_SetMinimumEpsilonStep(GPPropagatorInField *This,
						    G4double newEpsMin );

FQUALIFIER 
G4double  GPPropagatorInField_GetMaximumEpsilonStep( GPPropagatorInField *This );

FQUALIFIER 
void      GPPropagatorInField_SetMaximumEpsilonStep(GPPropagatorInField *This,
						    G4double newEpsMax );

FQUALIFIER
void GPPropagatorInField_SetLargestAcceptableStep(GPPropagatorInField *This,
						  G4double newBigDist );

FQUALIFIER
G4double  GPPropagatorInField_GetLargestAcceptableStep( GPPropagatorInField *This );

FQUALIFIER
GPFieldManager*  GPPropagatorInField_GetCurrentFieldManager( GPPropagatorInField *This );

FQUALIFIER
void GPPropagatorInField_SetThresholdNoZeroStep( GPPropagatorInField *This,
						 G4int noAct,
						 G4int noHarsh,
						 G4int noAbandon );

FQUALIFIER
G4int GPPropagatorInField_GetThresholdNoZeroSteps( GPPropagatorInField *This,
						   G4int i );

FQUALIFIER 
G4double  GPPropagatorInField_GetZeroStepThreshold( GPPropagatorInField *This );

FQUALIFIER 
void      GPPropagatorInField_SetZeroStepThreshold(GPPropagatorInField *This,
						   G4double newLength );

FQUALIFIER
void GPPropagatorInField_SetDetectorFieldManager(GPPropagatorInField *This,
						 GPFieldManager* newDetectorFieldManager);

FQUALIFIER
void GPPropagatorInField_SetUseSafetyForOptimization(GPPropagatorInField *This,
						     G4bool value );

FQUALIFIER 
G4bool 
GPPropagatorInField_GetUseSafetyForOptimization( GPPropagatorInField *This ); 

FQUALIFIER 
void GPPropagatorInField_SetNavigatorForPropagating(GPPropagatorInField *This, 
						    GPNavigator *SimpleOrMultiNavigator );

FQUALIFIER
GPNavigator* GPPropagatorInField_GetNavigatorForPropagating( GPPropagatorInField *This );

FQUALIFIER 
void GPPropagatorInField_SetIntersectionLocator(GPPropagatorInField *This,
						GPMultiLevelLocator *pIntLoc );

FQUALIFIER
GPMultiLevelLocator* GPPropagatorInField_GetIntersectionLocator( GPPropagatorInField *This );

FQUALIFIER
G4bool GPPropagatorInField_IntersectChord( GPPropagatorInField *This,
					   GPThreeVector StartPointA, 
					   GPThreeVector EndPointB,
					   G4double      &NewSafety,
					   G4double      &LinearStepLength,
					   GPThreeVector &IntersectionPoint);

}

#endif 
