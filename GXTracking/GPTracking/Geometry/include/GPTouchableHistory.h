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
// $Id: G4TouchableHistory.hh,v 1.11 2009-11-06 11:10:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4TouchableHistory
//
// Class description:
//
// Object representing a touchable detector element, and its history in the
// geometrical hierarchy, including its net resultant local->global transform.

// History:
// - Created. Paul Kent, August 1996
// ----------------------------------------------------------------------
#ifndef GPTOUCHABLEHISTORY_HH
#define GPTOUCHABLEHISTORY_HH

#include "GPNavigationHistory.h"
#include "GPThreeVector.h"
#include "GPRotationMatrix.h"

struct GPTouchableHistory
{
  GPRotationMatrix frot;
  GPThreeVector ftlate;
  GPNavigationHistory fhistory;
};

extern "C" {
FQUALIFIER
void GPTouchableHistory_Constructor(GPTouchableHistory *This);

FQUALIFIER
GPThreeVector GPTouchableHistory_GetTranslation(GPTouchableHistory *This,
						G4int depth);

FQUALIFIER
GPRotationMatrix GPTouchableHistory_GetRotation(GPTouchableHistory *This,
						G4int depth);

FQUALIFIER
void  GPTouchableHistory_UpdateYourself( GPTouchableHistory *This,
					 GPVPhysicalVolume*  pPhysVol,
					 GPNavigationHistory* pHistory ); 

FQUALIFIER
G4int GPTouchableHistory_CalculateHistoryIndex( GPTouchableHistory *This,
						G4int stackDepth );

FQUALIFIER
GPVPhysicalVolume* GPTouchableHistory_GetVolume( GPTouchableHistory *This,
						 G4int depth );

FQUALIFIER
GPVSolid* GPTouchableHistory_GetSolid( GPTouchableHistory *This,
				       G4int depth );

FQUALIFIER
G4int GPTouchableHistory_GetReplicaNumber( GPTouchableHistory *This,
					   G4int depth );

FQUALIFIER
G4int GPTouchableHistory_GetHistoryDepth( GPTouchableHistory *This);

FQUALIFIER
G4int GPTouchableHistory_MoveUpHistory( GPTouchableHistory *This,
					G4int num_levels );

FQUALIFIER
GPNavigationHistory* GPTouchableHistory_GetHistory( GPTouchableHistory *This );

FQUALIFIER
G4int GPTouchableHistory_GetCopyNumber(GPTouchableHistory *This, 
				       G4int depth);

}

#endif
