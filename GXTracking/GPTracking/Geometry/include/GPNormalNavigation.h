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
// $Id: G4NormalNavigation.hh,v 1.6 2010-11-04 08:57:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4NormalNavigation
//
// Class description:
//
// Utility for navigation in volumes containing only G4PVPlacement
// daughter volumes.

// History:
// - Created. Paul Kent, Aug 96
// --------------------------------------------------------------------
#ifndef GPNORMALNAVIGATION_HH
#define GPNORMALNAVIGATION_HH

//#include <iomanip>

#include "GPNavigationHistory.h"
//#include "G4NavigationLogger.hh"
#include "GPVPhysicalVolume.h"
#include "GPLogicalVolume.h"
#include "GPVSolid.h"
#include "GPThreeVector.h"
//#include "G4AuxiliaryNavServices.hh"    // Needed for inline methods

struct GPNormalNavigation
{
    G4bool fCheck; 
  //    G4NavigationLogger* fLogger;
};

extern "C" {

FQUALIFIER
void GPNormalNavigation_GPNormalNavigation(GPNormalNavigation *This);

FQUALIFIER
G4double
GPNormalNavigation_ComputeStep( GPNormalNavigation *This,
				GPThreeVector localPoint,
                                GPThreeVector localDirection,
                                G4double currentProposedStepLength,
				G4double *newSafety,
				GPNavigationHistory *history,
				G4bool *validExitNormal,
				GPThreeVector *exitNormal,
				G4bool *exiting,
				G4bool *entering,
				GEOMETRYLOC GPVPhysicalVolume *(*pBlockedPhysical),
				G4int *blockedReplicaNo);

FQUALIFIER
G4double GPNormalNavigation_ComputeSafety(GPNormalNavigation *This,
					  GPThreeVector localPoint,
					  GPNavigationHistory *history,
					  const G4double);

FQUALIFIER
G4bool GPNormalNavigation_LevelLocate( GPNormalNavigation *This,
				       GPNavigationHistory *history,
				       GEOMETRYLOC GPVPhysicalVolume* blockedVol,
				       GPThreeVector  globalPoint,
				       GPThreeVector* globalDirection,
				       G4bool  pLocatedOnEdge, 
				       GPThreeVector* localPoint );

FQUALIFIER
G4int GPNormalNavigation_GetVerboseLevel(GPNormalNavigation *This);

FQUALIFIER
void GPNormalNavigation_CheckMode(GPNormalNavigation *This, 
				  G4bool mode);

}
#endif
