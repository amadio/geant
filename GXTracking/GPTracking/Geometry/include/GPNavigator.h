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
// $Id: G4Navigator.hh,v 1.34 2010-12-15 13:46:39 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4Navigator
//
// Class description:
//
// A class for use by the tracking management, able to obtain/calculate
// dynamic tracking time information such as the distance to the next volume,
// or to find the physical volume containing a given point in the world
// reference system. The navigator maintains a transformation history and
// other information to optimise the tracking time performance.
// History:
// - Created.                                  Paul Kent,     Jul 95/96
// - Zero step protections                     J.A. / G.C.,   Nov  2004
// - Added check mode                          G. Cosmo,      Mar  2004
// - Made Navigator Abstract                   G. Cosmo,      Nov  2003
// *********************************************************************

#ifndef GPNAVIGATOR_HH
#define GPNAVIGATOR_HH 1

#include "GPTypeDef.h"
#include "GPGeomdefs.h"

#include "GPThreeVector.h"
#include "GPAffineTransform.h"
#include "GPRotationMatrix.h"

#include "GPVSolid.h"
#include "GPLogicalVolume.h"             // Used in inline methods

//#include "G4GRSVolume.hh"                 //    "         "
//#include "G4GRSSolid.hh"                  //    "         "
//#include "G4TouchableHandle.hh"           //    "         "
//#include "G4TouchableHistoryHandle.hh"

//#include "G4ParameterisedNavigation.hh"
//#include "G4ReplicaNavigation.hh"
//#include "G4RegularNavigation.hh"

//#include <iostream>

#include "GPVPhysicalVolume.h"
#include "GPNavigationHistory.h"
#include "GPNormalNavigation.h"
#include "GPVoxelNavigation.h"
#include "GPSaveNavigatorState.h"
#include "GPTouchableHistory.h"

struct GPNavigator
{
  G4double kCarTolerance;
  // Geometrical tolerance for surface thickness of shapes.
  //
  // BEGIN State information
  //
  GPNavigationHistory fHistory;
  // Transformation and history of the current path
  // through the geometrical hierarchy.

  G4bool fEnteredDaughter;
  // A memory of whether in this Step a daughter volume is entered 
  // (set in Compute & Locate).
  //  After Compute: it expects to enter a daughter
  //  After Locate:  it has entered a daughter
  
  G4bool fExitedMother;
  // A similar memory whether the Step exited current "mother" volume
  // completely, not entering daughter.

  G4bool fWasLimitedByGeometry;
  // Set true if last Step was limited by geometry.
  
  GPThreeVector fStepEndPoint;
  // Endpoint of last ComputeStep 
  // can be used for optimisation (e.g. when computing safety).
  GPThreeVector fLastStepEndPointLocal; 
  // Position of the end-point of the last call to ComputeStep 
  // in last Local coordinates.
  
  G4int  fVerbose;
  // Verbose(ness) level  [if > 0, printout can occur].
  
  G4bool fActive;
  // States if the navigator is activated or not.
  
  G4bool fLastTriedStepComputation; 
  // Whether ComputeStep was called since the last call to a Locate method
  // Uses: - distinguish parts of state which differ before/after calls
  //         to ComputeStep or one of the Locate methods;
  //       - avoid two consecutive calls to compute-step (illegal).
  
  G4bool fEntering,fExiting;
  // Entering/Exiting volumes blocking/setup
  // o If exiting
  //      volume ptr & replica number (set & used by Locate..())
  //      used for blocking on redescent of geometry
  // o If entering
  //      volume ptr & replica number (set by ComputeStep(),used by
  //      Locate..()) of volume for `automatic' entry
  
  GEOMETRYLOC GPVPhysicalVolume *fBlockedPhysicalVolume;
  G4int fBlockedReplicaNo;

  GPThreeVector fLastLocatedPointLocal;
  // Position of the last located point relative to its containing volume.
  G4bool fLocatedOutsideWorld;
  // Whether the last call to Locate methods left the world
  
  G4bool fValidExitNormal;    // Set true if have leaving volume normal
  GPThreeVector fExitNormal;  // Leaving volume normal, in the
                              // volume containing the exited
                              // volume's coordinate system
  GPThreeVector fGrandMotherExitNormal;  // Leaving volume normal, in its 
                                         // own coordinate system
  
  // Count zero steps - as one or two can occur due to changing momentum at
  //                    a boundary or at an edge common between volumes
  //                  - several are likely a problem in the geometry
  //                    description or in the navigation
  //
  G4bool fLastStepWasZero;
  // Whether the last ComputeStep moved Zero. Used to check for edges.
  
  G4bool fLocatedOnEdge;       
  // Whether the Navigator has detected an edge
  G4int fNumberZeroSteps;
  // Number of preceding moves that were Zero. Reset to 0 after finite step
  G4int fActionThreshold_NoZeroSteps;  
  // After this many failed/zero steps, act (push etc) 
  G4int fAbandonThreshold_NoZeroSteps; 
  // After this many failed/zero steps, abandon track
  
  GPThreeVector  fPreviousSftOrigin;
  G4double       fPreviousSafety; 
  // Memory of last safety origin & value. Used in ComputeStep to ensure
  // that origin of current Step is in the same volume as the point of the
  // last relocation
  
  //
  // END State information
  //
  
  // Save key state information (NOT the navigation history stack)
  //

  // Tracking Invariants
  //
  GEOMETRYLOC GPVPhysicalVolume  *fTopPhysical;
    // A link to the topmost physical volume in the detector.
    // Must be positioned at the origin and unrotated.

  GPSaveNavigatorState fSaveState;

  // Utility information
  //
  G4bool fCheck;
    // Check-mode flag  [if true, more strict checks are performed].
  G4bool fPushed, fWarnPush;
    // Push flags  [if true, means a stuck particle has been pushed].

  // Helpers/Utility classes
  //
  GPNormalNavigation  fnormalNav;
  GPVoxelNavigation fvoxelNav;
  //  G4ParameterisedNavigation fparamNav;
  //  G4ReplicaNavigation freplicaNav;
  //  G4RegularNavigation fregularNav;
};

extern "C" {

FQUALIFIER 
void GPNavigator_Constructor(GPNavigator *This);

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume*
GPNavigator_ResetHierarchyAndLocate(GPNavigator *This,
                                    GPThreeVector p,
                                    GPThreeVector *direction,
                                    GPTouchableHistory* h);

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* 
GPNavigator_LocateGlobalPointAndSetup( GPNavigator *This,
                                       GPThreeVector globalPoint,
                                       GPThreeVector* pGlobalDirection,
                                       const G4bool relativeSearch,
                                       const G4bool ignoreDirection );

FQUALIFIER
void
GPNavigator_LocateGlobalPointWithinVolume(GPNavigator *This,
                                          GPThreeVector pGlobalpoint);

FQUALIFIER
void GPNavigator_SetSavedState(GPNavigator *This);

FQUALIFIER
void GPNavigator_RestoreSavedState(GPNavigator *This);

FQUALIFIER
G4double GPNavigator_ComputeStep( GPNavigator *This,
                                  GPThreeVector pGlobalpoint,
                                  GPThreeVector pDirection,
                                  const G4double pCurrentProposedStepLength,
                                  G4double *pNewSafety);

FQUALIFIER
G4double GPNavigator_CheckNextStep(GPNavigator *This, 
                                   GPThreeVector pGlobalpoint,
                                   GPThreeVector pDirection,
                                   const G4double pCurrentProposedStepLength,
                                   G4double *pNewSafety);

FQUALIFIER
void GPNavigator_ResetState(GPNavigator *This);

FQUALIFIER
void GPNavigator_SetupHierarchy(GPNavigator *This);

FQUALIFIER
GPThreeVector GPNavigator_GetLocalExitNormal(GPNavigator *This,
                                             G4bool* valid );

FQUALIFIER
GPAffineTransform
GPNavigator_GetMotherToDaughterTransform( GPNavigator *This,
                                          GPVPhysicalVolume *pEnteringPhysVol,
                                          G4int   enteringReplicaNo,
                                          EVolume enteringVolumeType );

FQUALIFIER
GPThreeVector 
GPNavigator_GetLocalExitNormalAndCheck(GPNavigator *This,
				       GPThreeVector ExpectedBoundaryPointGlobal,
                                       G4bool*        pValid);
FQUALIFIER
GPThreeVector 
GPNavigator_GetGlobalExitNormal(GPNavigator *This,
                                GPThreeVector IntersectPointGlobal,
                                G4bool*        pValidNormal);

FQUALIFIER
G4double GPNavigator_ComputeSafety(GPNavigator *This, 
                                   GPThreeVector pGlobalpoint,
                                   const G4double pMaxLength,
                                   const G4bool keepState);
FQUALIFIER
GPThreeVector GPNavigator_GetCurrentLocalCoordinate(GPNavigator *This);

FQUALIFIER
GPThreeVector GPNavigator_ComputeLocalAxis(GPNavigator *This, 
                                           GPThreeVector pVec);

FQUALIFIER
GPThreeVector
GPNavigator_ComputeLocalPoint(GPNavigator *This,
                              GPThreeVector pGlobalPoint);

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* GPNavigator_GetWorldVolume(GPNavigator *This);

FQUALIFIER
void GPNavigator_SetWorldVolume(GPNavigator *This,
                                GEOMETRYLOC GPVPhysicalVolume* pWorld);

FQUALIFIER
void GPNavigator_SetGeometricallyLimitedStep(GPNavigator *This);

FQUALIFIER
void GPNavigator_ResetStackAndState(GPNavigator *This);

FQUALIFIER
EVolume GPNavigator_VolumeType(GPNavigator *This,
                               GEOMETRYLOC GPVPhysicalVolume *pVol);

FQUALIFIER
EVolume GPNavigator_CharacteriseDaughters(GPNavigator *This,
                                          GEOMETRYLOC GPLogicalVolume *pLog);

FQUALIFIER
GPAffineTransform GPNavigator_GetGlobalToLocalTransform(GPNavigator *This);

FQUALIFIER
GPAffineTransform  GPNavigator_GetLocalToGlobalTransform(GPNavigator *This);

FQUALIFIER
GPThreeVector GPNavigator_NetTranslation(GPNavigator *This);

FQUALIFIER
GPRotationMatrix GPNavigator_NetRotation(GPNavigator *This);

FQUALIFIER
G4bool GPNavigator_IsActive(GPNavigator *This);

FQUALIFIER
void GPNavigator_Activate(GPNavigator *This,
                          G4bool flag);

FQUALIFIER
G4bool GPNavigator_EnteredDaughterVolume(GPNavigator *This);

FQUALIFIER
G4bool GPNavigator_ExitedMotherVolume(GPNavigator *This);

FQUALIFIER
G4bool GPNavigator_IsCheckModeActive(GPNavigator *This);

FQUALIFIER
void GPNavigator_SetPushVerbosity(GPNavigator *This,
                                  G4bool mode);
FQUALIFIER 
G4int GPNavigator_SeverityOfZeroStepping(GPNavigator *This,
                                         G4int* noZeroSteps );
FQUALIFIER 
void GPNavigator_EnableBestSafety(GPNavigator *This,
                                  G4bool value );

}

#endif
