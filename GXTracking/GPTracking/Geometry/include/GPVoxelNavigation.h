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
// $Id: G4VoxelNavigation.hh,v 1.8 2010-11-04 12:13:30 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4VoxelNavigation
//
// Class description:
//
// Utility for navigation in volumes containing only G4PVPlacement
// daughter volumes for which voxels have been constructed.

// History:
// - Created. Paul Kent, Aug 96
// --------------------------------------------------------------------
#ifndef GPVOXELNAVIGATION_HH
#define GPVOXELNAVIGATION_HH

#include "GPGeomdefs.h"
#include "GPNavigationHistory.h"
//#include "G4NavigationLogger.hh"
#include "GPAffineTransform.h"
#include "GPVPhysicalVolume.h"
#include "GPLogicalVolume.h"
#include "GPVSolid.h"
#include "GPThreeVector.h"

//#include "G4BlockingList.hh"

//class G4VoxelSafety; 

// Required for inline implementation
//
//#include "G4AuxiliaryNavServices.hh"

// Required for voxel handling & voxel stack
//
//#include <vector>
#include "GPSmartVoxelProxy.h"
#include "GPSmartVoxelNode.h"
#include "GPSmartVoxelHeader.h"

//#ifndef __CUDA_ARCH__
//#define USE_BLIST
//#include <string.h>
//#endif

#define K_MAX_VOXEL_STACK_DEPTH 4 // mind alignment

struct GPVoxelNavigation
{
    //    G4BlockingList fBList;
    // Blocked volumes

    //
    //  BEGIN Voxel Stack information
    //

    G4int fVoxelDepth;
      // Note: fVoxelDepth==0+ => fVoxelAxisStack(0+) contains axes of voxel
      //       fVoxelDepth==-1 -> not in voxel

    EAxis fVoxelAxisStack[K_MAX_VOXEL_STACK_DEPTH];
      // Voxel axes

    G4int fVoxelNoSlicesStack[K_MAX_VOXEL_STACK_DEPTH];
      // No slices per voxel at each level

    G4double fVoxelSliceWidthStack[K_MAX_VOXEL_STACK_DEPTH]; 
      // Width of voxels at each level 

    G4int fVoxelNodeNoStack[K_MAX_VOXEL_STACK_DEPTH];    
      // Node no point is inside at each level 

    GEOMETRYLOC GPSmartVoxelHeader* fVoxelHeaderStack[K_MAX_VOXEL_STACK_DEPTH];
      // Voxel headers at each level

    GEOMETRYLOC GPSmartVoxelNode* fVoxelNode;
      // Node containing last located point

    //
    //  END Voxel Stack information
    //

    //    G4VoxelSafety  *fpVoxelSafety;
    //  Helper object for Voxel Safety

    //    G4bool fCheck;
    //    G4bool fBestSafety; 

    //    G4NavigationLogger* fLogger;
    // Verbosity logger

  //#ifdef USE_BLIST
  //    char *fBlist;
  //    int fBlistSz;
  //#endif

};

extern "C" {

FQUALIFIER
void GPVoxelNavigation_Constructor( GPVoxelNavigation *This );

FQUALIFIER
G4double GPVoxelNavigation_ComputeVoxelSafety(GPVoxelNavigation *This,
					      GPThreeVector localPoint);

FQUALIFIER
G4bool GPVoxelNavigation_LevelLocate(GPVoxelNavigation *This,
				     GPNavigationHistory* history,
				     GEOMETRYLOC const GPVPhysicalVolume* blockedVol,
				     GPThreeVector globalPoint,
				     GPThreeVector* globalDirection,
				     const G4bool pLocatedOnEdge, 
				     GPThreeVector *localPoint );

FQUALIFIER
G4double GPVoxelNavigation_ComputeStep(GPVoxelNavigation *This,
				       GPThreeVector localPoint,
				       GPThreeVector localDirection,
				       const G4double currentProposedStepLength,
				       G4double *newSafety,
				       GPNavigationHistory *history,
				       G4bool *validExitNormal,
				       GPThreeVector *exitNormal,
				       G4bool *exiting,
				       G4bool *entering,
				       GEOMETRYLOC GPVPhysicalVolume *(*pBlockedPhysical) );

FQUALIFIER
GEOMETRYLOC GPSmartVoxelNode* GPVoxelNavigation_VoxelLocate(GPVoxelNavigation *This,
						GEOMETRYLOC GPSmartVoxelHeader* pHead,
						GPThreeVector localPoint);

FQUALIFIER
G4bool GPVoxelNavigation_LocateNextVoxel(GPVoxelNavigation *This,
					 GPThreeVector localPoint,
					 GPThreeVector localDirection,
					 const G4double currentStep );


FQUALIFIER
G4double GPVoxelNavigation_ComputeSafety(GPVoxelNavigation *This,
					 GPThreeVector localPoint,
					 GPNavigationHistory *history ,
					 const G4double       maxLength);

}

#endif
