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
// $Id: G4NavigationHistory.hh,v 1.19 2010-12-15 17:05:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4NavigationHistory
//
// Class description:
//
// Responsible for maintenance of the history of the path taken through
// the geometrical hierarchy. Principally a utility class for use by the
// G4Navigator.

// History:
//
// 25.07.96 P.Kent Initial version. Services derived from
//                 requirements of G4Navigator.
// ----------------------------------------------------------------------
#ifndef GPNAVIGATIONHISTORY_HH
#define GPNAVIGATIONHISTORY_HH

//#include <assert.h>
#include "GPGeomdefs.h"

#include "GPAffineTransform.h"
#include "GPVPhysicalVolume.h"
#include "GPNavigationLevel.h"
//#include "G4EnhancedVecAllocator.hh"

//#include <vector>
//#include <iostream>

struct GPNavigationHistory
{
  //  std::vector<G4NavigationLevel> fNavHistory;
  GPNavigationLevel fNavHistory[kHistoryMax];
  G4int fStackDepth;
    // Depth of stack: effectively depth in geometrical tree

};

extern "C" {

FQUALIFIER
void GPNavigationHistory_Constructor(GPNavigationHistory *This);

FQUALIFIER
void GPNavigationHistory_Reset(GPNavigationHistory *This);

FQUALIFIER
void GPNavigationHistory_Clear(GPNavigationHistory *This);

FQUALIFIER
void GPNavigationHistory_SetFirstEntry(GPNavigationHistory *This,
				       GPVPhysicalVolume* pVol);

FQUALIFIER
const GPAffineTransform* 
GPNavigationHistory_GetPtrTopTransform(GPNavigationHistory *This);

FQUALIFIER
GPAffineTransform 
GPNavigationHistory_GetTopTransform(GPNavigationHistory *This);

FQUALIFIER
G4int GPNavigationHistory_GetTopReplicaNo(GPNavigationHistory *This);

FQUALIFIER
EVolume GPNavigationHistory_GetTopVolumeType(GPNavigationHistory *This);

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* GPNavigationHistory_GetTopVolume(GPNavigationHistory *This);

FQUALIFIER
G4int GPNavigationHistory_GetDepth(GPNavigationHistory *This);

FQUALIFIER
GPAffineTransform GPNavigationHistory_GetTransform(GPNavigationHistory *This,
						   G4int n);
FQUALIFIER
G4int GPNavigationHistory_GetReplicaNo(GPNavigationHistory *This,
				       G4int n);

FQUALIFIER
EVolume GPNavigationHistory_GetVolumeType(GPNavigationHistory *This,
					  G4int n);

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* GPNavigationHistory_GetVolume(GPNavigationHistory *This,
						 G4int n);

FQUALIFIER
G4int GPNavigationHistory_GetMaxDepth(GPNavigationHistory *This);

FQUALIFIER
void GPNavigationHistory_BackLevel(GPNavigationHistory *This);

FQUALIFIER
void GPNavigationHistory_BackLevel2(GPNavigationHistory *This,
				    G4int n);

FQUALIFIER
void GPNavigationHistory_NewLevel( GPNavigationHistory *This,
				   GEOMETRYLOC GPVPhysicalVolume *pNewMother,
				   EVolume vType,
				   G4int nReplica );
}

#endif
