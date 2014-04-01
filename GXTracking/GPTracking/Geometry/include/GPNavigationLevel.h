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
// $Id: G4NavigationLevel.hh,v 1.18 2009-08-04 08:27:20 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4NavigationLevel
//
// Class description:
//
// Maintains one level of the geometrical hierarchy.
// A utility class for use by G4NavigationHistory.

// History:
//
// 30.09.97 J.Apostolakis Initial version. Services derived from
//                        requirements of touchables & G4NavigatorHistory.
// ----------------------------------------------------------------------
#ifndef GPNAVIGATIONLEVEL_HH
#define GPNAVIGATIONLEVEL_HH

#include "GPTypeDef.h"

#include "GPAffineTransform.h"
#include "GPVPhysicalVolume.h"

//#include "GPNavigationLevelRep.h"
//#include "G4Allocator.hh"

struct GPNavigationLevel
{
  //   GPNavigationLevelRep*  fLevelRep;
   GPAffineTransform  sTransform;
     // Compounded global->local transformation (takes a point in the 
     // global reference system to the system of the volume at this level)

   GEOMETRYLOC GPVPhysicalVolume* sPhysicalVolumePtr;
     // Physical volume ptrs, for this level's volume

   G4int              sReplicaNo;
   EVolume            sVolumeType;
     // Volume `type' 

   G4int              fCountRef; 

};

extern "C" {

FQUALIFIER
void GPNavigationLevel_Constructor( GPNavigationLevel *This,
				    GEOMETRYLOC GPVPhysicalVolume* pPhysVol,
				    GPAffineTransform  afTransform,
				    EVolume            volTp,
				    G4int              repNo );

FQUALIFIER
void GPNavigationLevel_Constructor2( GPNavigationLevel *This,
				     GEOMETRYLOC GPVPhysicalVolume *pPhysVol,
				     GPAffineTransform levelAbove,
				     GPAffineTransform relativeCurrent,
				     EVolume            volTp,
				     G4int              repNo );

FQUALIFIER
void GPNavigationLevel_GPNavigationLevel0(GPNavigationLevel *This);

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* GPNavigationLevel_GetPhysicalVolume(GPNavigationLevel *This);

FQUALIFIER
GPAffineTransform GPNavigationLevel_GetTransform(GPNavigationLevel *This);

FQUALIFIER
GPAffineTransform* GPNavigationLevel_GetPtrTransform(GPNavigationLevel *This); 

FQUALIFIER
EVolume GPNavigationLevel_GetVolumeType(GPNavigationLevel *This);

FQUALIFIER
G4int GPNavigationLevel_GetReplicaNo(GPNavigationLevel *This);

FQUALIFIER
void GPNavigationLevel_AddAReference(GPNavigationLevel *This); 

FQUALIFIER
G4bool GPNavigationLevel_RemoveAReference(GPNavigationLevel *This);
 
}

#endif
