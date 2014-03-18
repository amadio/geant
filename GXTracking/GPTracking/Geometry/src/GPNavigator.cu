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
// $Id: G4Navigator.cc,v 1.46 2010-11-15 14:03:27 gcosmo Exp $
// GEANT4 tag $ Name:  $
// 
// class G4Navigator Implementation
//
// Original author: Paul Kent, July 95/96
//
// --------------------------------------------------------------------

#include "GPNavigator.h"

//#include "G4ios.hh"
//#include <iomanip>
//#include "G4GeometryTolerance.hh"
//#include "GPVPhysicalVolume.h"

// ********************************************************************
// Constructor
// ********************************************************************
//
FQUALIFIER 
void GPNavigator_Constructor(GPNavigator *This)
{
  GPNavigationHistory_Constructor( &(This->fHistory) );

  This->fWasLimitedByGeometry=false;

  This->fVerbose=0;
  This->fTopPhysical=GEOMETRYNULL;
  This->fCheck=false;
  This->fPushed=false; 
  This->fWarnPush=true;

  This->fActive= false; 
  This->fLastTriedStepComputation= false;

  GPNavigator_ResetStackAndState(This);

  This->fActionThreshold_NoZeroSteps  = 10; 
  This->fAbandonThreshold_NoZeroSteps = 25; 

  This->kCarTolerance = 1E-9;
  //  GPNavigator_SetNormalNavigation(This->fregularNav, &fnormalNav );

  This->fStepEndPoint = GPThreeVector_create( kInfinity, kInfinity, kInfinity ); 
  This->fLastStepEndPointLocal = GPThreeVector_create( kInfinity, kInfinity, kInfinity ); 

  // this->SetVerboseLevel(3);
  // this->CheckMode(true);
}

// ********************************************************************
// Destructor
// ********************************************************************
//
//GPNavigator_~G4Navigator()
//{;}

// ********************************************************************
// ResetHierarchyAndLocate
// ********************************************************************
//
FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume*
GPNavigator_ResetHierarchyAndLocate(GPNavigator *This,
				    GPThreeVector p,
				    GPThreeVector *direction,
				    GPTouchableHistory* h
)
{
  GPNavigator_ResetState(This);
  This->fHistory = *(GPTouchableHistory_GetHistory(h));
  GPNavigator_SetupHierarchy(This);
  This->fLastTriedStepComputation= false;  // Redundant, but best
  return GPNavigator_LocateGlobalPointAndSetup(This, p, direction, true, false);
}

// ********************************************************************
// LocateGlobalPointAndSetup
//
// Locate the point in the hierarchy return 0 if outside
// The direction is required 
//    - if on an edge shared by more than two surfaces 
//      (to resolve likely looping in tracking)
//    - at initial location of a particle
//      (to resolve potential ambiguity at boundary)
// 
// Flags on exit: (comments to be completed)
// fEntering         - True if entering `daughter' volume (or replica)
//                     whether daughter of last mother directly 
//                     or daughter of that volume's ancestor.
// ********************************************************************
//
FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* 
GPNavigator_LocateGlobalPointAndSetup( GPNavigator *This,
				       GPThreeVector globalPoint,
				       GPThreeVector* pGlobalDirection,
				       const G4bool relativeSearch,
				       const G4bool ignoreDirection )
{
  G4bool notKnownContained=true, noResult;
  GPVPhysicalVolume *targetPhysical;
  GPLogicalVolume *targetLogical;
  GPVSolid *targetSolid=GEOMETRYNULL;

  GPThreeVector localPoint = GPThreeVector_create(0,0,0);
  GPThreeVector globalDirection = GPThreeVector_create(0,0,0);
  EInside insideCode;
  
  G4bool considerDirection = (!ignoreDirection) || This->fLocatedOnEdge;
  This->fLastTriedStepComputation= false;   

  if( considerDirection && pGlobalDirection != 0 )
  {
    globalDirection= *pGlobalDirection;
  }

  if ( !relativeSearch )
  {
    GPNavigator_ResetStackAndState(This);
  }
  else
  {
    if ( This->fWasLimitedByGeometry )
    {
      This->fWasLimitedByGeometry = false;
      This->fEnteredDaughter = This->fEntering;   // Remember
      This->fExitedMother = This->fExiting;       // Remember
      if ( This->fExiting )
      {
        if ( GPNavigationHistory_GetDepth(&(This->fHistory)) )
        {
          This->fBlockedPhysicalVolume = GPNavigationHistory_GetTopVolume(&(This->fHistory));
          This->fBlockedReplicaNo = GPNavigationHistory_GetTopReplicaNo(&(This->fHistory));
          GPNavigationHistory_BackLevel(&(This->fHistory));
        }
        else
        {
          This->fLastLocatedPointLocal = localPoint;
          This->fLocatedOutsideWorld = true;
          return GEOMETRYNULL;           // Have exited world volume
        }
        // A fix for the case where a volume is "entered" at an edge
        // and a coincident surface exists outside it.
        //  - This stops it from exiting further volumes and cycling
        //  - However ReplicaNavigator treats this case itself
        //
        if ( This->fLocatedOnEdge && 
	 (GPNavigator_VolumeType(This,This->fBlockedPhysicalVolume)!=kReplica ))
        { 
          This->fExiting= false;
        }
      }
      else
        if ( This->fEntering )
        {
          switch (GPNavigator_VolumeType(This,This->fBlockedPhysicalVolume))
          {
            case kNormal:
              GPNavigationHistory_NewLevel(&(This->fHistory),
		                  This->fBlockedPhysicalVolume, kNormal,0);
	      //GPVPhysicalVolume_GetCopyNo(This->fBlockedPhysicalVolume));
              break;
	      //@@@G4FWP
	      /*
            case kReplica:
              freplicaNav.ComputeTransformation(fBlockedReplicaNo,
                                                fBlockedPhysicalVolume);
              fHistory.NewLevel(fBlockedPhysicalVolume, kReplica,
                                fBlockedReplicaNo);
              fBlockedPhysicalVolume->SetCopyNo(fBlockedReplicaNo);
              break;
            case kParameterised:
              if( fBlockedPhysicalVolume->GetRegularStructureId() == 0 )
              {
                GPVSolid *pSolid;
                G4VPVParameterisation *pParam;
                G4TouchableHistory parentTouchable( fHistory );
                pParam = fBlockedPhysicalVolume->GetParameterisation();
                pSolid = pParam->ComputeSolid(fBlockedReplicaNo,
                                              fBlockedPhysicalVolume);
                pSolid->ComputeDimensions(pParam, fBlockedReplicaNo,
                                          fBlockedPhysicalVolume);
                pParam->ComputeTransformation(fBlockedReplicaNo,
                                              fBlockedPhysicalVolume);
                fHistory.NewLevel(fBlockedPhysicalVolume, kParameterised,
                                  fBlockedReplicaNo);
                fBlockedPhysicalVolume->SetCopyNo(fBlockedReplicaNo);
                //
                // Set the correct solid and material in Logical Volume
                //
                GPLogicalVolume *pLogical;
                pLogical = fBlockedPhysicalVolume->GetLogicalVolume();
                pLogical->SetSolid( pSolid );
                pLogical->UpdateMaterial(pParam ->
                  ComputeMaterial(fBlockedReplicaNo,
                                  fBlockedPhysicalVolume, 
                                  &parentTouchable));
              }
              break;
	      */
          }
          This->fEntering = false;
          This->fBlockedPhysicalVolume = GEOMETRYNULL;
          GPAffineTransform aT = GPNavigationHistory_GetTopTransform( &(This->fHistory) );
          localPoint =  GPAffineTransform_TransformPoint(&aT,globalPoint);
          notKnownContained = false;
        }
    }
    else
    {
      This->fBlockedPhysicalVolume = GEOMETRYNULL;
      This->fEntering = false;
      This->fEnteredDaughter = false;  // Full Step was not taken, did not enter
      This->fExiting = false;
      This->fExitedMother = false;     // Full Step was not taken, did not exit
    }
  }

  //
  // Search from top of history up through geometry until
  // containing volume found:
  // If on 
  // o OUTSIDE - Back up level, not/no longer exiting volumes
  // o SURFACE and EXITING - Back up level, setting new blocking no.s
  // else
  // o containing volume found
  //

  while (notKnownContained)
  {
    if ( GPNavigationHistory_GetTopVolumeType(&(This->fHistory))!=kReplica )
    {
      //      targetSolid = fHistory.GetTopVolume()->GetLogicalVolume()->GetSolid();

      GPVPhysicalVolume*  testv = GPNavigationHistory_GetTopVolume(&(This->fHistory));
      GPLogicalVolume*  testLog = GPVPhysicalVolume_GetLogicalVolume(testv);

      targetSolid = GPLogicalVolume_GetSolid(
		    GPVPhysicalVolume_GetLogicalVolume(
		    GPNavigationHistory_GetTopVolume(&(This->fHistory))));

      GPAffineTransform aT = GPNavigationHistory_GetTopTransform( &(This->fHistory) );

      localPoint = GPAffineTransform_TransformPoint(&aT,globalPoint);
      insideCode = GPVSolid_Inside(targetSolid,localPoint);
    }
    else
    {
      //@@@G4FWP no replica
      //      insideCode = GPNavigator_BackLocate(This->freplicaNav, &(This->fHistory), globalPoint, localPoint,
      //                                          This->fExiting, notKnownContained);
      //@@@G4FWP no replica

      // !CARE! if notKnownContained returns false then the point is within
      // the containing placement volume of the replica(s). If insidecode
      // will result in the history being backed up one level, then the
      // local point returned is the point in the system of this new level
    }

    if ( insideCode==kOutside )
    {

      if ( GPNavigationHistory_GetDepth(&(This->fHistory)) )
      {
        This->fBlockedPhysicalVolume =  GPNavigationHistory_GetTopVolume(&(This->fHistory));
        This->fBlockedReplicaNo =  GPNavigationHistory_GetTopReplicaNo(&(This->fHistory));
        GPNavigationHistory_BackLevel(&(This->fHistory));
        This->fExiting = false;
      }
      else
      {
        This->fLastLocatedPointLocal = localPoint;
        This->fLocatedOutsideWorld = true;
        return GEOMETRYNULL;         // Have exited world volume
      }
    }
    else

      if ( insideCode==kSurface )
      {
        G4bool isExiting = This->fExiting;
        if( (!(This->fExiting))&&considerDirection )
        {
          // Figure out whether we are exiting this level's volume
          // by using the direction
          //
          G4bool directionExiting = false;

	  GPAffineTransform aT = GPNavigationHistory_GetTopTransform( &(This->fHistory) );
          GPThreeVector localDirection = GPAffineTransform_TransformAxis(&aT,globalDirection);
          if ( GPNavigationHistory_GetTopVolumeType(&(This->fHistory))!=kReplica )
          {
            GPThreeVector normal = GPVSolid_SurfaceNormal(targetSolid,
							  localPoint);
            directionExiting = GPThreeVector_dot(normal,localDirection) > 0.0;
            isExiting = isExiting || directionExiting;
          }
        }
        if( isExiting )
        {
          if ( GPNavigationHistory_GetDepth(&(This->fHistory)) )
          {
            This->fBlockedPhysicalVolume =  GPNavigationHistory_GetTopVolume(&(This->fHistory));
            This->fBlockedReplicaNo =   GPNavigationHistory_GetTopReplicaNo(&(This->fHistory));
            GPNavigationHistory_BackLevel(&(This->fHistory));
            //
            // Still on surface but exited volume not necessarily convex
            //
            This->fValidExitNormal = false;
          } 
          else
          {
            This->fLastLocatedPointLocal = localPoint;
            This->fLocatedOutsideWorld = true;
            return GEOMETRYNULL;          // Have exited world volume
          }
        }
        else
        {
          notKnownContained=false;
        }
      }
      else
      {
        notKnownContained=false;
      }
  }  // END while (notKnownContained)

  //
  // Search downwards until deepest containing volume found,
  // blocking fBlockedPhysicalVolume/BlockedReplicaNum
  //
  // 3 Cases:
  //
  // o Parameterised daughters
  //   =>Must be one G4PVParameterised daughter & voxels
  // o Positioned daughters & voxels
  // o Positioned daughters & no voxels

  noResult = true;  // noResult should be renamed to 
                    // something like enteredLevel, as that is its meaning.
  do
  {
    // Determine `type' of current mother volume
    //
    targetPhysical = GPNavigationHistory_GetTopVolume(&(This->fHistory));
    if (!targetPhysical) { break; }
    targetLogical = GPVPhysicalVolume_GetLogicalVolume(targetPhysical);
    switch( GPNavigator_CharacteriseDaughters(This,targetLogical) )
    {
      case kNormal:

        if ( GPLogicalVolume_GetVoxelHeader(targetLogical) != GEOMETRYNULL )  // use optimised navigation
        {
          noResult = GPVoxelNavigation_LevelLocate(&(This->fvoxelNav),
						   &(This->fHistory),
						   This->fBlockedPhysicalVolume,
						   globalPoint,
						   pGlobalDirection,
						   considerDirection,
						   &localPoint);
        }
        else                       // do not use optimised navigation
	{

	  noResult = GPNormalNavigation_LevelLocate(&(This->fnormalNav),
						    &(This->fHistory),
						    This->fBlockedPhysicalVolume,
						    globalPoint,
						    pGlobalDirection,
						    considerDirection,
						    &localPoint);
	}
        break;
	/*
      case kReplica:
        noResult = freplicaNav.LevelLocate(fHistory,
                                           fBlockedPhysicalVolume,
                                           fBlockedReplicaNo,
                                           globalPoint,
                                           pGlobalDirection,
                                           considerDirection,
                                           localPoint);
        break;
      case kParameterised:
        if( GetDaughtersRegularStructureId(targetLogical) != 1 )
        {
          noResult = fparamNav.LevelLocate(fHistory,
                                           fBlockedPhysicalVolume,
                                           fBlockedReplicaNo,
                                           globalPoint,
                                           pGlobalDirection,
                                           considerDirection,
                                           localPoint);
        }
        else  // Regular structure
        {
          noResult = fregularNav.LevelLocate(fHistory,
                                             fBlockedPhysicalVolume,
                                             fBlockedReplicaNo,
                                             globalPoint,
                                             pGlobalDirection,
                                             considerDirection,
                                             localPoint);
        }
        break;
	*/

    }

    // LevelLocate returns true if it finds a daughter volume 
    // in which globalPoint is inside (or on the surface).

    if ( noResult )
    {
      // Entering a daughter after ascending
      //
      // The blocked volume is no longer valid - it was for another level
      //
      This->fBlockedPhysicalVolume = GEOMETRYNULL;
      This->fBlockedReplicaNo = -1;

      // fEntering should be false -- else blockedVolume is assumed good.
      // fEnteredDaughter is used for ExitNormal
      //
      This->fEntering = false;
      This->fEnteredDaughter = true;
    }
  } while (noResult);

  This->fLastLocatedPointLocal = localPoint;

  This->fLocatedOutsideWorld= false;

  return targetPhysical;
}

// ********************************************************************
// LocateGlobalPointWithinVolume
//
// -> the state information of this Navigator and its subNavigators
//    is updated in order to start the next step at pGlobalpoint
// -> no check is performed whether pGlobalpoint is inside the 
//    original volume (this must be the case).
//
// Note: a direction could be added to the arguments, to aid in future
//       optional checking (via the old code below, flagged by OLD_LOCATE). 
//       [ This would be done only in verbose mode ]
// ********************************************************************
//
FQUALIFIER
void
GPNavigator_LocateGlobalPointWithinVolume(GPNavigator *This,
					  GPThreeVector pGlobalpoint)
{  
  This->fLastLocatedPointLocal = GPNavigator_ComputeLocalPoint(This,
							       pGlobalpoint);
  This->fLastTriedStepComputation= false;

   // For the case of Voxel (or Parameterised) volume the respective 
   // Navigator must be messaged to update its voxel information etc

   // Update the state of the Sub Navigators 
   // - in particular any voxel information they store/cache
   //
   GPVPhysicalVolume*  motherPhysical = GPNavigationHistory_GetTopVolume(&(This->fHistory));
   GPLogicalVolume*    motherLogical  = GPVPhysicalVolume_GetLogicalVolume(motherPhysical);
   //@@@G4FWP no voxel
   //   G4SmartVoxelHeader* pVoxelHeader = GPLogicalVolume_GetVoxelHeader(motherLogical);

   if ( GPNavigationHistory_GetTopVolumeType(&(This->fHistory))!=kReplica )
   {
     switch( GPNavigator_CharacteriseDaughters(This,motherLogical) )
     {
       //@@@G4FWP No VolelNavigation
       /*
       case kNormal:
         if ( pVoxelHeader )
         {
           GPVoxelNavigation_VoxelLocate(&(This->fVoxelNav), pVoxelHeader, This->fLastLocatedPointLocal );
         }
         break;
       case kParameterised:
         if( GetDaughtersRegularStructureId(motherLogical) != 1 )
         {
           // Resets state & returns voxel node
           //
           fparamNav.ParamVoxelLocate( pVoxelHeader, fLastLocatedPointLocal );
         }
         break;
       case kReplica:
	 ;
         break;
       */
     }
   }

   // Reset the state variables 
   //   - which would have been affected
   //     by the 'equivalent' call to LocateGlobalPointAndSetup
   //   - who's values have been invalidated by the 'move'.
   //
   This->fBlockedPhysicalVolume = GEOMETRYNULL; 
   This->fBlockedReplicaNo = -1;
   This->fEntering = false;
   This->fEnteredDaughter = false;  // Boundary not encountered, did not enter
   This->fExiting = false;
   This->fExitedMother = false;     // Boundary not encountered, did not exit
}

// ********************************************************************
// SetSavedState
//
// Save the state, in case this is a parasitic call
// Save fValidExitNormal, fExitNormal, fExiting, fEntering, 
//      fBlockedPhysicalVolume, fBlockedReplicaNo, fLastStepWasZero; 
// ********************************************************************
//
FQUALIFIER
void GPNavigator_SetSavedState(GPNavigator *This)
{
  // fSaveExitNormal = fExitNormal; 
  (This->fSaveState).sExitNormal = This->fExitNormal;
  (This->fSaveState).sValidExitNormal = This->fValidExitNormal;
  (This->fSaveState).sExiting = This->fExiting;
  (This->fSaveState).sEntering = This->fEntering;

  (This->fSaveState).spBlockedPhysicalVolume = This->fBlockedPhysicalVolume;
  (This->fSaveState).sBlockedReplicaNo = This->fBlockedReplicaNo; 
  
  (This->fSaveState).sLastStepWasZero = This->fLastStepWasZero; 
}

// ********************************************************************
// RestoreSavedState
//
// Restore the state (in Compute Step), in case this is a parasitic call
// ********************************************************************
//
FQUALIFIER
void GPNavigator_RestoreSavedState(GPNavigator *This)
{
  This->fExitNormal = (This->fSaveState).sExitNormal;
  This->fValidExitNormal = (This->fSaveState).sValidExitNormal;
  This->fExiting = (This->fSaveState).sExiting;
  This->fEntering = (This->fSaveState).sEntering;

  This->fBlockedPhysicalVolume = (This->fSaveState).spBlockedPhysicalVolume;
  This->fBlockedReplicaNo = (This->fSaveState).sBlockedReplicaNo; 

  This->fLastStepWasZero = (This->fSaveState).sLastStepWasZero; 
}

// ********************************************************************
// ComputeStep
//
// Computes the next geometric Step: intersections with current
// mother and `daughter' volumes.
//
// NOTE:
//
// Flags on entry:
// --------------
// fValidExitNormal  - Normal of exited volume is valid (convex, not a 
//                     coincident boundary)
// fExitNormal       - Surface normal of exited volume
// fExiting          - True if have exited solid
//
// fBlockedPhysicalVolume - Ptr to exited volume (or 0)
// fBlockedReplicaNo - Replication no of exited volume
// fLastStepWasZero  - True if last Step size was zero.
//
// Flags on exit:
// -------------
// fValidExitNormal  - True if surface normal of exited volume is valid
// fExitNormal       - Surface normal of exited volume rotated to mothers
//                    reference system
// fExiting          - True if exiting mother
// fEntering         - True if entering `daughter' volume (or replica)
// fBlockedPhysicalVolume - Ptr to candidate (entered) volume
// fBlockedReplicaNo - Replication no of candidate (entered) volume
// fLastStepWasZero  - True if this Step size was zero.
// ********************************************************************
//
FQUALIFIER
G4double GPNavigator_ComputeStep( GPNavigator *This,
				  GPThreeVector pGlobalpoint,
				  GPThreeVector pDirection,
				  const G4double pCurrentProposedStepLength,
				  G4double *pNewSafety)
{
  GPThreeVector localDirection = GPNavigator_ComputeLocalAxis(This,pDirection);
  G4double Step = kInfinity;
  GEOMETRYLOC GPVPhysicalVolume  *motherPhysical = GPNavigationHistory_GetTopVolume(&(This->fHistory));
  GEOMETRYLOC GPLogicalVolume *motherLogical = GPVPhysicalVolume_GetLogicalVolume(motherPhysical);

  //  static G4int sNavCScalls=0;
  //  sNavCScalls++;

  This->fLastTriedStepComputation= true; 

  GPThreeVector newLocalPoint = GPNavigator_ComputeLocalPoint(This,pGlobalpoint);
  if( GPThreeVector_nequal(newLocalPoint,This->fLastLocatedPointLocal) )
  {
    // Check whether the relocation is within safety
    //
    GPThreeVector oldLocalPoint = This->fLastLocatedPointLocal;
    G4double moveLenSq = GPThreeVector_mag2(GPThreeVector_sub(newLocalPoint,oldLocalPoint));

    if ( moveLenSq >= This->kCarTolerance*This->kCarTolerance )
    {
      // Relocate the point within the same volume
      //
      GPNavigator_LocateGlobalPointWithinVolume(This, pGlobalpoint );
      This->fLastTriedStepComputation= true;     // Ensure that this is set again !!
    }
  }

  if ( GPNavigationHistory_GetTopVolumeType(&(This->fHistory))!=kReplica ) {
    switch( GPNavigator_CharacteriseDaughters(This,motherLogical) ) {
    case kNormal:
      
      if ( GPLogicalVolume_GetVoxelHeader(motherLogical) != GEOMETRYNULL ) {
	Step = GPVoxelNavigation_ComputeStep(&(This->fvoxelNav),
					     This->fLastLocatedPointLocal,
					     localDirection,
					     pCurrentProposedStepLength,
					     pNewSafety,
					     &(This->fHistory),
					     &(This->fValidExitNormal),
					     &(This->fExitNormal),
					     &(This->fExiting),
					     &(This->fEntering),
					     &(This->fBlockedPhysicalVolume));
      }
      else {
	//@@@G4FWP
	// if( GPVPhysicalVolume_GetRegularStructureId(motherPhysical) == 0 ) {
	// Returns non-zero code in case the underlying volume structure 
        //  is regular, voxel-like. 
        if(1) {
	  //@@@G4FWP - check this line is necessary
	  // GPNavigator_LocateGlobalPointAndSetup(This, pGlobalpoint, NULL, 
	  //					   false, true );
	  
	  Step = GPNormalNavigation_ComputeStep(&(This->fnormalNav),
						This->fLastLocatedPointLocal,
						localDirection,
						pCurrentProposedStepLength,
						pNewSafety,
						&(This->fHistory),
						&(This->fValidExitNormal),
						&(This->fExitNormal),
						&(This->fExiting),
						&(This->fEntering),
						&(This->fBlockedPhysicalVolume),
						&(This->fBlockedReplicaNo));
	}
	else { // Regular (non-voxelised) structure

	  GPNavigator_LocateGlobalPointAndSetup(This, pGlobalpoint, &pDirection, true, true );
	  This->fLastTriedStepComputation= true;     // Ensure that this is set again !!
	  //
	  // if physical process limits the step, the voxel will not be the
	  // one given by ComputeStepSkippingEqualMaterials() and the local
	  // point will be wrongly calculated.
	  
	  // There is a problem: when msc limits the step and the point is
	  // assigned wrongly to phantom in previous step (while it is out
	  // of the container volume). Then LocateGlobalPointAndSetup() has
	  // reset the history topvolume to world.
	  //
	  //@@@G4FWP never happen
	  //if(GPVPhysicalVolume_GetRegularStructureId(GPNavigationHistory_GetTopVolume(&(This->fHistory)) ) == 0 )
	  if( 0 ) { 
	    //              G4Exception("GPNavigator_ComputeStep()",
	    //                          "GeomNav1001", JustWarning,
	    //                "Point is relocated in voxels, while it should be outside!");
	    Step = GPNormalNavigation_ComputeStep(&(This->fnormalNav),
						  This->fLastLocatedPointLocal,
						  localDirection,
						  pCurrentProposedStepLength,
						  pNewSafety,
						  &(This->fHistory),
						  &(This->fValidExitNormal),
						  &(This->fExitNormal),
						  &(This->fExiting),
						  &(This->fEntering),
						  &(This->fBlockedPhysicalVolume),
						  &(This->fBlockedReplicaNo));
	  }
	  /*
	    else
	    {
	    Step = fregularNav.
	    ComputeStepSkippingEqualMaterials(fLastLocatedPointLocal,
	    localDirection,
	    pCurrentProposedStepLength,
	    pNewSafety,
	    fHistory,
	    fValidExitNormal,
	    fExitNormal,
	    fExiting,
	    fEntering,
	    &fBlockedPhysicalVolume,
	    fBlockedReplicaNo,
	    motherPhysical);
	    }
	  */
	}
      }
      break;
      /*
	case kParameterised:
	if( GetDaughtersRegularStructureId(motherLogical) != 1 )
	{
	Step = fparamNav.ComputeStep(fLastLocatedPointLocal,
	localDirection,
	pCurrentProposedStepLength,
	pNewSafety,
	fHistory,
	fValidExitNormal,
	fExitNormal,
	fExiting,
	fEntering,
	&fBlockedPhysicalVolume,
	fBlockedReplicaNo);
	}
	else  // Regular structure
	{
	Step = fregularNav.ComputeStep(fLastLocatedPointLocal,
	localDirection,
	pCurrentProposedStepLength,
	pNewSafety,
	fHistory,
	fValidExitNormal,
	fExitNormal,
	fExiting,
	fEntering,
	&fBlockedPhysicalVolume,
	fBlockedReplicaNo);
	}
	break;
	case kReplica:
	//        G4Exception("GPNavigator_ComputeStep()", "GeomNav0001",
	//        FatalException, "Not applicable for replicated volumes.");
	break;
      */
    }
  }  
  /*

  }

  else
  {
    // In the case of a replica, it must handle the exiting
    // edge/corner problem by itself
    //
    G4bool exitingReplica = fExitedMother;
    Step = freplicaNav.ComputeStep(pGlobalpoint,
                                   pDirection,
                                   fLastLocatedPointLocal,
                                   localDirection,
                                   pCurrentProposedStepLength,
                                   pNewSafety,
                                   fHistory,
                                   fValidExitNormal,
                                   fExitNormal,
                                   exitingReplica,
                                   fEntering,
                                   &fBlockedPhysicalVolume,
                                   fBlockedReplicaNo);
    fExiting= exitingReplica;                          // still ok to set it ??
  }
  */

  // Remember last safety origin & value.
  //
  This->fPreviousSftOrigin = pGlobalpoint;
  This->fPreviousSafety = *pNewSafety; 

  // Count zero steps - one can occur due to changing momentum at a boundary
  //                  - one, two (or a few) can occur at common edges between
  //                    volumes
  //                  - more than two is likely a problem in the geometry
  //                    description or the Navigation 

  // Rule of thumb: likely at an Edge if two consecutive steps are zero,
  //                because at least two candidate volumes must have been
  //                checked
  //
  This->fLocatedOnEdge   = This->fLastStepWasZero && (Step==0.0);
  This->fLastStepWasZero = (Step==0.0);
  if (This->fPushed)  This->fPushed = This->fLastStepWasZero;

  // Handle large number of consecutive zero steps
  //
  if ( This->fLastStepWasZero )
  {
    This->fNumberZeroSteps++;
    if( This->fNumberZeroSteps > This->fActionThreshold_NoZeroSteps-1 )
    {
       // Act to recover this stuck track. Pushing it along direction
       //
       Step += 100*kCarTolerance;
       This->fPushed = true;
    }
    if( This->fNumberZeroSteps > This->fAbandonThreshold_NoZeroSteps-1 )
    {
      // Must kill this stuck track
      //
      //@@@G4FWP this is always false
      //      GPVPhysicalVolume_CheckOverlaps(1,5000, false);
    }
  }
  else
  {
    if (!This->fPushed)  This->fNumberZeroSteps = 0;
  }

  This->fEnteredDaughter = This->fEntering;   // I expect to enter a volume in this Step
  This->fExitedMother = This->fExiting;

  This->fStepEndPoint = GPThreeVector_add(pGlobalpoint,GPThreeVector_mult(pDirection,Step)); 
  This->fLastStepEndPointLocal = GPThreeVector_add(This->fLastLocatedPointLocal,
						   GPThreeVector_mult(localDirection,Step)); 

  if( This->fExiting )
  {
    if(This->fValidExitNormal)
    {
      // Convention: fExitNormal is in the 'grand-mother' coordinate system
      //
      This->fGrandMotherExitNormal= This->fExitNormal;
    }
    else
    {  
      // We must calculate the normal anyway (in order to have it if requested)
      //
      GPThreeVector finalLocalPoint = 
                  GPThreeVector_add(This->fLastLocatedPointLocal,
		  GPThreeVector_mult(localDirection,Step));
      // Now fGrandMotherExitNormal is in the 'grand-mother' coordinate system
      //
      This->fGrandMotherExitNormal = GPVSolid_SurfaceNormal(
				     GPLogicalVolume_GetSolid(motherLogical),
				     finalLocalPoint);

      const GPRotationMatrix* mRot = GPVPhysicalVolume_GetRotation(motherPhysical);
      if( mRot )
      { 
        GPRotationMatrix inv = GPRotationMatrix_inverse(mRot);
        This->fGrandMotherExitNormal = GPRotationMatrix_apply(&inv,This->fGrandMotherExitNormal);
	  //  *= (*mRot).inverse();
      }
      //  Do not set fValidExitNormal -- this signifies that the solid is convex!
    }
  }
  This->fStepEndPoint= GPThreeVector_add(pGlobalpoint,GPThreeVector_mult(pDirection,Step)); 

  if( (Step == pCurrentProposedStepLength) && (!(This->fExiting)) && (!(This->fEntering)) )
  {
    // This if Step is not really limited by the geometry.
    // The Navigator is obliged to return "infinity"
    //
    Step = kInfinity;
  }

  return Step;
}

// ********************************************************************
// CheckNextStep
//
// Compute the step without altering the navigator state
// ********************************************************************
//
FQUALIFIER
G4double GPNavigator_CheckNextStep(GPNavigator *This, 
				   GPThreeVector pGlobalpoint,
				   GPThreeVector pDirection,
				   const G4double pCurrentProposedStepLength,
				   G4double *pNewSafety)
{
  G4double step;

  // Save the state, for this parasitic call
  //
  GPNavigator_SetSavedState(This);

  step = GPNavigator_ComputeStep (This, 
				  pGlobalpoint, 
				  pDirection,
				  pCurrentProposedStepLength, 
				  pNewSafety ); 

  // If a parasitic call, then attempt to restore the key parts of the state
  //
  GPNavigator_RestoreSavedState(This); 

  return step; 
}

// ********************************************************************
// ResetState
//
// Resets stack and minimum of navigator state `machine'
// ********************************************************************
//
FQUALIFIER
void GPNavigator_ResetState(GPNavigator *This)
{
  This->fWasLimitedByGeometry  = false;
  This->fEntering              = false;
  This->fExiting               = false;
  This->fLocatedOnEdge         = false;
  This->fLastStepWasZero       = false;
  This->fEnteredDaughter       = false;
  This->fExitedMother          = false;
  This->fPushed                = false;

  This->fValidExitNormal       = false;
  This->fExitNormal            = GPThreeVector_create(0,0,0);

  This->fPreviousSftOrigin     = GPThreeVector_create(0,0,0);
  This->fPreviousSafety        = 0.0; 

  This->fNumberZeroSteps       = 0;
    
  This->fBlockedPhysicalVolume = GEOMETRYNULL;
  This->fBlockedReplicaNo      = -1;

  This->fLastLocatedPointLocal = GPThreeVector_create( kInfinity, -kInfinity, 0.0 ); 
  This->fLocatedOutsideWorld   = false;
}

// ********************************************************************
// SetupHierarchy
//
// Renavigates & resets hierarchy described by current history
// o Reset volumes
// o Recompute transforms and/or solids of replicated/parameterised volumes
// ********************************************************************
//
FQUALIFIER
void GPNavigator_SetupHierarchy(GPNavigator *This)
{
  G4int i;
  const G4int cdepth = GPNavigationHistory_GetDepth(&(This->fHistory));
  //  GPVPhysicalVolume *current;
  //  GPVSolid *pSolid;
  //  GPVPVParameterisation *pParam;

  for ( i=1; i<=cdepth; i++ )
  {
    //    current = GPNavigationHistory_GetVolume(&(This->fHistory), i);
    switch (  GPNavigationHistory_GetVolumeType(&(This->fHistory), i) )
    {
      case kNormal:
        break;
	/*
      case kReplica:
        freplicaNav.ComputeTransformation(fHistory.GetReplicaNo(i), current);
        break;
      case kParameterised:
        G4int replicaNo;
        pParam = current->GetParameterisation();
        replicaNo = fHistory.GetReplicaNo(i);
        pSolid = pParam->ComputeSolid(replicaNo, current);

        // Set up dimensions & transform in solid/physical volume
        //
        pSolid->ComputeDimensions(pParam, replicaNo, current);
        pParam->ComputeTransformation(replicaNo, current);

        G4TouchableHistory touchable( fHistory );
        touchable.MoveUpHistory();  // move up to the parent level
      
        // Set up the correct solid and material in Logical Volume
        //
        GPLogicalVolume *pLogical = current->GetLogicalVolume();
        pLogical->SetSolid( pSolid );
        pLogical->UpdateMaterial( pParam ->
          ComputeMaterial(replicaNo, current, &touchable) );
        break;
	*/
    }
  }
}

// ********************************************************************
// GetLocalExitNormal
//
// Obtains the Normal vector to a surface (in local coordinates)
// pointing out of previous volume and into current volume
// ********************************************************************
//
FQUALIFIER
GPThreeVector GPNavigator_GetLocalExitNormal(GPNavigator *This,
					     G4bool* valid )
{
  GPThreeVector    ExitNormal = GPThreeVector_create(0.,0.,0.);
  GPVSolid        *currentSolid=0;
  GPLogicalVolume *candidateLogical;

  if ( This->fLastTriedStepComputation ) 
  {
    // use fLastLocatedPointLocal
    // and next candidate volume
    GPThreeVector nextSolidExitNormal = GPThreeVector_create(0.,0.,0.);

    if( This->fEntering && ((This->fBlockedPhysicalVolume)!=0) ) 
    { 
      candidateLogical= GPVPhysicalVolume_GetLogicalVolume(This->fBlockedPhysicalVolume);
      if( candidateLogical ) 
      {
        // fLastStepEndPointLocal is in the coordinates of the mother
        // we need it in the daughter's coordinate system.

        if( GPNavigator_CharacteriseDaughters(This,candidateLogical) != kReplica )
        {
          // First transform fLastLocatedPointLocal to the new daughter
          // coordinates
          GPAffineTransform MotherToDaughterTransform=
            GPNavigator_GetMotherToDaughterTransform(This, 
						     This->fBlockedPhysicalVolume, 
						     This->fBlockedReplicaNo,
						     GPNavigator_VolumeType(This,This->fBlockedPhysicalVolume) ); 
          GPThreeVector daughterPointOwnLocal= 
            GPAffineTransform_TransformPoint(&MotherToDaughterTransform,This->fLastStepEndPointLocal ); 

          // OK if it is a parameterised volume
          //
          EInside  inSideIt; 
          G4bool   onSurface;
          G4double safety= -1.0; 

          currentSolid= GPLogicalVolume_GetSolid(candidateLogical); 
          inSideIt  =   GPVSolid_Inside(currentSolid,daughterPointOwnLocal); 

          onSurface =  (inSideIt == kSurface); 
          if( ! onSurface ) 
          {
            if( inSideIt == kOutside )
            { 
              safety = (GPVSolid_DistanceToIn(currentSolid,daughterPointOwnLocal)); 
              onSurface = safety < 100.0 * kCarTolerance; 
            }
            else if (inSideIt == kInside ) 
            {
              safety = (GPVSolid_DistanceToOut(currentSolid,daughterPointOwnLocal)); 
              onSurface = safety < 100.0 * kCarTolerance; 
            }
          }

          if( onSurface ) 
          {
            nextSolidExitNormal =
              GPVSolid_SurfaceNormal(currentSolid,daughterPointOwnLocal); 
 
            // Entering the solid ==> opposite
            //
            ExitNormal = GPThreeVector_mult(nextSolidExitNormal,-1.0);
          }
          else
          {
          }
          *valid = onSurface;   //   was =true;
        }
        else
        {
          *valid = false;  // TODO: Need Separate code for replica!!!!
        }
      }
    }
    else if ( This->fExiting ) 
    {
        ExitNormal = This->fGrandMotherExitNormal;
        *valid = true;
    }
    else  // ie  ( fBlockedPhysicalVolume == 0 )
    {
      *valid = false;
    }
  }
  else 
  {
    if ( GPNavigator_EnteredDaughterVolume(This) )
    {
      ExitNormal= GPThreeVector_mult(
                  GPVSolid_SurfaceNormal(
		  GPLogicalVolume_GetSolid(
                  GPVPhysicalVolume_GetLogicalVolume(
                  GPNavigationHistory_GetTopVolume(&(This->fHistory)))),
                  This->fLastLocatedPointLocal),-1.0);

      *valid = true;
    }
    else
    {
      if( This->fExitedMother )
      {
        ExitNormal = This->fGrandMotherExitNormal;
        *valid = true;
      }
      else  // We are not at a boundary. ExitNormal remains (0,0,0)
      {
        *valid = false;
      }
    }
  }
  return ExitNormal;
}

// ********************************************************************
// GetMotherToDaughterTransform
//
// Obtains the mother to daughter affine transformation
// ********************************************************************
//
FQUALIFIER
GPAffineTransform
GPNavigator_GetMotherToDaughterTransform( GPNavigator *This,
					  GEOMETRYLOC GPVPhysicalVolume *pEnteringPhysVol,   // not Const 
					  G4int   enteringReplicaNo,
					  EVolume enteringVolumeType ) 
{
  switch (enteringVolumeType)
  {
    case kNormal:  // Nothing is needed to prepare the transformation
      break;       // It is stored already in the physical volume (placement)
      /*
    case kReplica: // Sets the transform in the Replica - tbc
      //      G4Exception("GPNavigator_GetMotherToDaughterTransform()",
      //                  "GeomNav0001", FatalException,
      //                  "Method NOT Implemented yet for replica volumes.");
      break;
    case kParameterised:
      if( pEnteringPhysVol->GetRegularStructureId() == 0 )
      {
        G4VPVParameterisation *pParam =
          pEnteringPhysVol->GetParameterisation();
        GPVSolid* pSolid =
          pParam->ComputeSolid(enteringReplicaNo, pEnteringPhysVol);
        pSolid->ComputeDimensions(pParam, enteringReplicaNo, pEnteringPhysVol);

        // Sets the transform in the Parameterisation
        //
        pParam->ComputeTransformation(enteringReplicaNo, pEnteringPhysVol);

        // Set the correct solid and material in Logical Volume
        //
        GPLogicalVolume* pLogical = pEnteringPhysVol->GetLogicalVolume();
        pLogical->SetSolid( pSolid );
      }
      break;
      */
  }
  GPAffineTransform  aT;
  GPAffineTransform_Constructor3(&aT,
				 GPVPhysicalVolume_GetRotation(pEnteringPhysVol), 
				 GPVPhysicalVolume_GetTranslation(pEnteringPhysVol));

  return GPAffineTransform_Invert(&aT); 
}

// ********************************************************************
// GetLocalExitNormalAndCheck
//
// Obtains the Normal vector to a surface (in local coordinates)
// pointing out of previous volume and into current volume, and
// checks the current point against expected 'local' value.
// ********************************************************************
//
FQUALIFIER
GPThreeVector 
GPNavigator_GetLocalExitNormalAndCheck(GPNavigator *This,
				       GPThreeVector ExpectedBoundaryPointGlobal,
				       G4bool*        pValid)
{
  GPThreeVector ExpectedBoundaryPointLocal = GPThreeVector_create(0,0,0);

  // Check Current point against expected 'local' value
  //
  if ( This->fLastTriedStepComputation ) 
  {
     GPAffineTransform GlobalToLocal= GPNavigator_GetGlobalToLocalTransform(This); 
     ExpectedBoundaryPointLocal = 
       GPAffineTransform_TransformPoint(&GlobalToLocal, ExpectedBoundaryPointGlobal ); 
  }

  return GPNavigator_GetLocalExitNormal(This, pValid); 
}

// ********************************************************************
// GetGlobalExitNormal
//
// Obtains the Normal vector to a surface (in global coordinates)
// pointing out of previous volume and into current volume
// ********************************************************************
//
FQUALIFIER
GPThreeVector 
GPNavigator_GetGlobalExitNormal(GPNavigator *This,
				GPThreeVector IntersectPointGlobal,
				G4bool*        pValidNormal)
{
  G4bool         validNormal;
  GPThreeVector  localNormal, globalNormal;

  localNormal = GPNavigator_GetLocalExitNormalAndCheck(This, IntersectPointGlobal, &validNormal);
  *pValidNormal = validNormal; 
  GPAffineTransform localToGlobal = GPNavigator_GetLocalToGlobalTransform(This); 
  globalNormal = GPAffineTransform_TransformAxis(&localToGlobal, localNormal );
  
  return globalNormal;
}

// ********************************************************************
// ComputeSafety
//
// It assumes that it will be 
//  i) called at the Point in the same volume as the EndPoint of the
//     ComputeStep.
// ii) after (or at the end of) ComputeStep OR after the relocation.
// ********************************************************************
//
FQUALIFIER
G4double GPNavigator_ComputeSafety(GPNavigator *This, 
				   GPThreeVector pGlobalpoint,
				   const G4double pMaxLength,
				   const G4bool keepState)
{
  G4double newSafety = 0.0;

  if (keepState)  { GPNavigator_SetSavedState(This); }
  //  fLastTriedStepComputation= true;   -- this method is NOT computing the Step size

  G4double distEndpointSq =  GPThreeVector_mag2(GPThreeVector_sub(pGlobalpoint,This->fStepEndPoint)); 
  G4bool   stayedOnEndpoint  = distEndpointSq < kCarTolerance*kCarTolerance; 
  G4bool   endpointOnSurface = This->fEnteredDaughter || This->fExitedMother;

  if( !(endpointOnSurface && stayedOnEndpoint) )
  {
    // Pseudo-relocate to this point (updates voxel information only)
    //
    GPNavigator_LocateGlobalPointWithinVolume(This, pGlobalpoint );
      // --->> Danger: Side effects on sub-navigator voxel information <<---
      //       Could be replaced again by 'granular' calls to sub-navigator
      //       locates (similar side-effects, but faster.  
      //       Solutions:
      //        1) Re-locate (to where?)
      //        2) Insure that the methods using (G4ComputeStep?)
      //           does a relocation (if information is disturbed only ?)

    GPVPhysicalVolume *motherPhysical = GPNavigationHistory_GetTopVolume(&(This->fHistory));
    GPLogicalVolume *motherLogical = GPVPhysicalVolume_GetLogicalVolume(motherPhysical);
    GPSmartVoxelHeader* pVoxelHeader = GPLogicalVolume_GetVoxelHeader(motherLogical);
    GPThreeVector localPoint = GPNavigator_ComputeLocalPoint(This,pGlobalpoint);

    if ( GPNavigationHistory_GetTopVolumeType(&(This->fHistory))!=kReplica )
    {
      switch(GPNavigator_CharacteriseDaughters(This,motherLogical))
      {
        case kNormal:
          if ( pVoxelHeader != GEOMETRYNULL )
          {
            newSafety= GPVoxelNavigation_ComputeSafety(&(This->fvoxelNav),localPoint,&(This->fHistory),pMaxLength);
          }
          else
          {
	    newSafety= GPNormalNavigation_ComputeSafety(&(This->fnormalNav),localPoint,&(This->fHistory),pMaxLength);
	  }
          break;
	  /*
        case kParameterised:
          if( GetDaughtersRegularStructureId(motherLogical) != 1 )
          {
            newSafety = fparamNav.ComputeSafety(localPoint,fHistory,pMaxLength);
          }
          else  // Regular structure
          {
            newSafety = fregularNav.ComputeSafety(localPoint,fHistory,pMaxLength);
          }
          break;
        case kReplica:
	  // G4Exception("GPNavigator_ComputeSafety()", "NotApplicable",
	  // FatalException, "Not applicable for replicated volumes.");
          break;
	  */
      }
    }
    /*
    else
    {
      newSafety = freplicaNav.ComputeSafety(pGlobalpoint, localPoint,
                                            fHistory, pMaxLength);
    }
    */
  }
  else // if( endpointOnSurface && stayedOnEndpoint )
  {
    newSafety = 0.0; 
  }

  // Remember last safety origin & value
  //
  This->fPreviousSftOrigin = pGlobalpoint;
  This->fPreviousSafety = newSafety; 

  if (keepState)  { GPNavigator_RestoreSavedState(This); }

  return newSafety;
}

// ********************************************************************
// CreateTouchableHistoryHandle
// ********************************************************************
//
//GPTouchableHistoryHandle GPNavigator_CreateTouchableHistoryHandle() 
//{
//  return G4TouchableHistoryHandle( CreateTouchableHistory() );
//}

// ********************************************************************
// PrintState
// ********************************************************************
//
FQUALIFIER
void  GPNavigator_PrintState()
{
  //do nothing
  ;
}

// ********************************************************************
// ComputeStepLog
// ********************************************************************
//
FQUALIFIER
void GPNavigator_ComputeStepLog(GPThreeVector pGlobalpoint,
				G4double moveLenSq)
{
  //do nothing
  ;

  //  The following checks only make sense if the move is larger
  //  than the tolerance.

  /*
  static const G4double fAccuracyForWarning   = kCarTolerance,
                        fAccuracyForException = 1000*kCarTolerance;

  GPThreeVector OriginalGlobalpoint = fHistory.GetTopTransform().Inverse().
                                      TransformPoint(fLastLocatedPointLocal); 

  G4double shiftOriginSafSq = (fPreviousSftOrigin-pGlobalpoint).mag2();

  // Check that the starting point of this step is 
  // within the isotropic safety sphere of the last point
  // to a accuracy/precision  given by fAccuracyForWarning.
  //   If so give warning.
  //   If it fails by more than fAccuracyForException exit with error.
  //
  if( shiftOriginSafSq >= sqr(fPreviousSafety) )
  {
    G4double shiftOrigin = std::sqrt(shiftOriginSafSq);
    G4double diffShiftSaf = shiftOrigin - fPreviousSafety;

    if( diffShiftSaf > fAccuracyForWarning )
    {
      G4int oldcoutPrec= G4cout.precision(8);
      G4int oldcerrPrec= G4cerr.precision(10);
      std::ostringstream message, suggestion;
      message << "Accuracy error or slightly inaccurate position shift."
              << G4endl
              << "     The Step's starting point has moved " 
              << std::sqrt(moveLenSq)/mm << " mm " << G4endl
              << "     since the last call to a Locate method." << G4endl
              << "     This has resulted in moving " 
              << shiftOrigin/mm << " mm " 
              << " from the last point at which the safety " 
              << "     was calculated " << G4endl
              << "     which is more than the computed safety= " 
              << fPreviousSafety/mm << " mm  at that point." << G4endl
              << "     This difference is " 
              << diffShiftSaf/mm << " mm." << G4endl
              << "     The tolerated accuracy is "
              << fAccuracyForException/mm << " mm.";

      suggestion << " ";
      static G4int warnNow = 0;
      if( ((++warnNow % 100) == 1) )
      {
        message << G4endl
               << "  This problem can be due to either " << G4endl
               << "    - a process that has proposed a displacement"
               << " larger than the current safety , or" << G4endl
               << "    - inaccuracy in the computation of the safety";
        suggestion << "We suggest that you " << G4endl
                   << "   - find i) what particle is being tracked, and "
                   << " ii) through what part of your geometry " << G4endl
                   << "      for example by re-running this event with "
                   << G4endl
                   << "         /tracking/verbose 1 "  << G4endl
                   << "    - check which processes you declare for"
                   << " this particle (and look at non-standard ones)"
                   << G4endl
                   << "   - in case, create a detailed logfile"
                   << " of this event using:" << G4endl
                   << "         /tracking/verbose 6 ";
      }
      G4Exception("GPNavigator_ComputeStep()",
                  "GeomNav1002", JustWarning,
                  message, G4String(suggestion.str()));
      G4cout.precision(oldcoutPrec);
      G4cerr.precision(oldcerrPrec);
    }
  }
  G4double safetyPlus = fPreviousSafety + fAccuracyForException;
  if ( shiftOriginSafSq > sqr(safetyPlus) )
  {
    std::ostringstream message;
    message << "May lead to a crash or unreliable results." << G4endl
            << "        Position has shifted considerably without"
            << " notifying the navigator !" << G4endl
            << "        Tolerated safety: " << safetyPlus << G4endl
            << "        Computed shift  : " << shiftOriginSafSq;
    G4Exception("GPNavigator_ComputeStep()", "GeomNav1002",
                JustWarning, message);
  }
  */
}

// ********************************************************************
// Operator <<
// ********************************************************************
//
//std::ostream& operator << (std::ostream &os,const G4Navigator &n)
//{
//  os << "Current History: " << G4endl << n.fHistory;
//  return os;
//}

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
// $Id: G4Navigator.icc,v 1.18 2010-12-15 13:46:39 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Navigator Inline implementation
//
// ********************************************************************

// ********************************************************************
// GetCurrentLocalCoordinate
//
// Returns the local coordinate of the current track
// ********************************************************************
//
FQUALIFIER
GPThreeVector GPNavigator_GetCurrentLocalCoordinate(GPNavigator *This)
{
  return This->fLastLocatedPointLocal;
}

// ********************************************************************
// ComputeLocalAxis
//
// Returns local direction of vector direction in world coord system
// ********************************************************************
//
FQUALIFIER
GPThreeVector GPNavigator_ComputeLocalAxis(GPNavigator *This, 
					   GPThreeVector pVec)
{
  GPAffineTransform aT =
                GPNavigationHistory_GetTopTransform(&(This->fHistory));

  return (  GPAffineTransform_IsRotated(&aT))
    ? GPAffineTransform_TransformAxis(&aT,pVec) : pVec ;
}

// ********************************************************************
// ComputeLocalPoint
//
// Returns local coordinates of a point in the world coord system
// ********************************************************************
//

FQUALIFIER
GPThreeVector
GPNavigator_ComputeLocalPoint(GPNavigator *This,
			      GPThreeVector pGlobalPoint)
{
  GPAffineTransform aT =
                GPNavigationHistory_GetTopTransform(&(This->fHistory));

  return ( GPAffineTransform_TransformPoint(&aT,pGlobalPoint) ) ;
}

// ********************************************************************
// GetWorldVolume
//
// Returns the current  world (`topmost') volume
// ********************************************************************
//
FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* GPNavigator_GetWorldVolume(GPNavigator *This)
{
  return This->fTopPhysical;
}

// ********************************************************************
// SetWorldVolume
//
// Sets the world (`topmost') volume
// ********************************************************************
//
FQUALIFIER
void GPNavigator_SetWorldVolume(GPNavigator *This,
				GEOMETRYLOC GPVPhysicalVolume* pWorld)
{
  //  if ( !(pWorld->GetTranslation()==GPThreeVector(0,0,0)) )
  //  {
  //    G4Exception ("GPNavigator_SetWorldVolume()", "GeomNav0002",
  //                 FatalException, "Volume must be centered on the origin.");
  //  }
  //  const GPRotationMatrix* rm = pWorld->GetRotation();
  //  if ( rm && (!rm->isIdentity()) )
  //  {
  //    G4Exception ("GPNavigator_SetWorldVolume()", "GeomNav0002",
  //                 FatalException, "Volume must not be rotated.");
  //  }

  This->fTopPhysical = pWorld;
  GPNavigationHistory_SetFirstEntry(&(This->fHistory),pWorld );

}

// ********************************************************************
// SetGeometrycallyLimitedStep
//
// Informs the navigator that the previous Step calculated
// by the geometry was taken in its entirety
// ********************************************************************
//
FQUALIFIER
void GPNavigator_SetGeometricallyLimitedStep(GPNavigator *This)
{
  This->fWasLimitedByGeometry=true;
}

// ********************************************************************
// ResetStackAndState
//
// Resets stack and minimum of navigator state `machine'
// ********************************************************************
//
FQUALIFIER
void GPNavigator_ResetStackAndState(GPNavigator *This)
{
  GPNavigationHistory_Reset(&(This->fHistory));
  GPNavigator_ResetState(This);
}

// ********************************************************************
// VolumeType
// ********************************************************************
//
FQUALIFIER
EVolume GPNavigator_VolumeType(GPNavigator *This,
			       GPVPhysicalVolume *pVol)
{
  EVolume type;

  //G4FWP
  /*
  //  EAxis axis;
  //  G4int nReplicas;
  //  G4double width,offset;
  //  G4bool consuming;
  if ( GPVPhysicalVolume_IsReplicated(pVol) )
  {
    GPVPhysicalVolume_GetReplicationData(pVol,axis,nReplicas,width,offset,consuming);
    type = (consuming) ? kReplica : kParameterised;
  }
  else
  {
    type = kNormal;
  }
  */  
  type = kNormal;
  return type;
}

// ********************************************************************
// CharacteriseDaughters
// ********************************************************************
//
FQUALIFIER
EVolume GPNavigator_CharacteriseDaughters(GPNavigator *This,
					  GEOMETRYLOC GPLogicalVolume *pLog)
{
  EVolume type;
  //  EAxis axis;
  //  G4int nReplicas;
  //  G4double width,offset;
  //  G4bool consuming;
  //  GPVPhysicalVolume pVol;

  if ( GPLogicalVolume_GetNoDaughters(pLog)==1 )
  {
    //    pVol = GPLogicalVolume_GetDaughter(pLog,0);

    /*
    if (GPVPhysicalVolume_IsReplicated(pVol))
    {
      GPVPhysicalVolume_GetReplicationData(pVol,axis,nReplicas,width,offset,consuming);
      type = (consuming) ? kReplica : kParameterised;
    }
    else
    {
      type = kNormal;
    }
    */
    type = kNormal;
  }
  else
  {
    type = kNormal;
  }
  return type;
}

// ********************************************************************
// GetDaughtersRegularStructureId
// ********************************************************************
//
/*
FQUALIFIER
G4int GPNavigator_GetDaughtersRegularStructureId(GPNavigator *This,
						 const GPLogicalVolume *pLog)
{
  G4int regId = 0;
  GPVPhysicalVolume *pVol;

  if ( GPLogicalVolume_GetNoDaughters(pLog)==1 )
  {
    pVol = GPLogicalVolume_GetDaughter(pLog,0);
    regId = GPVPhysicalVolume_GetRegularStructureId(pVol);
  }
  return regId;
}
*/

// ********************************************************************
// GetGlobalToLocalTransform
//
// Returns local to global transformation.
// I.e. transformation that will take point or axis in world coord system
// and return one in the local coord system
// ********************************************************************
//
FQUALIFIER
GPAffineTransform GPNavigator_GetGlobalToLocalTransform(GPNavigator *This)
{
  return GPNavigationHistory_GetTopTransform(&(This->fHistory));
}

// ********************************************************************
// GetLocalToGlobalTransform
//
// Returns global to local transformation 
// ********************************************************************
//
FQUALIFIER
GPAffineTransform  GPNavigator_GetLocalToGlobalTransform(GPNavigator *This) 
{
  GPAffineTransform  tempTransform;
  GPAffineTransform aT = GPNavigationHistory_GetTopTransform(&(This->fHistory));
  tempTransform = GPAffineTransform_Inverse(&aT); 
  return tempTransform;
}

// ********************************************************************
// NetTranslation
//
// Computes+returns the local->global translation of current volume
// ********************************************************************
//
FQUALIFIER
GPThreeVector GPNavigator_NetTranslation(GPNavigator *This)
{

  GPAffineTransform aT = GPNavigationHistory_GetTopTransform(&(This->fHistory));
  GPAffineTransform tf = GPAffineTransform_Inverse(&aT); 
  return GPAffineTransform_NetTranslation(&tf);
}

// ********************************************************************
// NetRotation
//
// Computes+returns the local->global rotation of current volume
// ********************************************************************
//
FQUALIFIER
GPRotationMatrix GPNavigator_NetRotation(GPNavigator *This)
{
 
  GPAffineTransform aT = GPNavigationHistory_GetTopTransform(&(This->fHistory));
  GPAffineTransform tf = GPAffineTransform_Inverse(&aT); 
  return GPAffineTransform_NetRotation(&tf);
}

// ********************************************************************
// CreateGRSVolume
//
// `Touchable' creation method: caller has deletion responsibility
// ********************************************************************
//
/*
G4GRSVolume* GPNavigator_CreateGRSVolume() const
{
  GPAffineTransform aT = GPNavigationHistory_GetTopTransform(&(This->fHistory));
  GPAffineTransform tf = GPAffineTransform_Inverse(&aT); 

  return new G4GRSVolume(fHistory.GetTopVolume(),
                         tf.NetRotation(),
                         tf.NetTranslation());
}
*/
// ********************************************************************
// CreateGRSSolid
//
// `Touchable' creation method: caller has deletion responsibility
// ********************************************************************
//
 /*
G4GRSSolid* GPNavigator_CreateGRSSolid() const
{
  GPAffineTransform tf(fHistory.GetTopTransform().Inverse());
  return new G4GRSSolid(fHistory.GetTopVolume()->GetLogicalVolume()->GetSolid(),
                        tf.NetRotation(),
                        tf.NetTranslation());
}
*/
// ********************************************************************
// CreateTouchableHistory
//
// `Touchable' creation method: caller has deletion responsibility
// ********************************************************************
//

/*
GPTouchableHistory* GPNavigator_CreateTouchableHistory(GPNavigator *This)
{
  //  return new GPTouchableHistory(This->fHistory);
  // use GPTouchableHistory_Constructor Directly
}
*/

// ********************************************************************
// CreateTouchableHistory(history)
//
// `Touchable' creation method: caller has deletion responsibility
// ********************************************************************
//
/*
GPTouchableHistory*
GPNavigator_CreateTouchableHistory(GPNavigator *This, 
GPNavigationHistory* history)
{
  //  return new G4TouchableHistory(*history);
  // using GPTouchableHistory_Constructor2 Directly
}
*/

// ********************************************************************
// LocateGlobalPointAndUpdateTouchableHandle
// ********************************************************************
//
  /*
void GPNavigator_LocateGlobalPointAndUpdateTouchableHandle(
                               const GPThreeVector&       position,
                               const GPThreeVector&       direction,
                                     G4TouchableHandle&   oldTouchableToUpdate,
                               const G4bool               RelativeSearch )
{
  GPVPhysicalVolume* pPhysVol;
  pPhysVol = LocateGlobalPointAndSetup( position,&direction,RelativeSearch );
  if( fEnteredDaughter || fExitedMother )
  {
     oldTouchableToUpdate = CreateTouchableHistory();
     if( pPhysVol == 0 )
     {
       // We want to ensure that the touchable is correct in this case.
       //  The method below should do this and recalculate a lot more ....
       //
       oldTouchableToUpdate->UpdateYourself( pPhysVol, &fHistory );
     }
  }
  return;
}
  */
// ********************************************************************
// LocateGlobalPointAndUpdateTouchable
//
// Use direction
// ********************************************************************
//
   /*
void GPNavigator_LocateGlobalPointAndUpdateTouchable(
                           const GPThreeVector&       position,
                           const GPThreeVector&       direction,
                                 G4VTouchable*        touchableToUpdate,
                           const G4bool               RelativeSearch  )
{
  GPVPhysicalVolume* pPhysVol;
  pPhysVol = LocateGlobalPointAndSetup( position, &direction, RelativeSearch);  
  touchableToUpdate->UpdateYourself( pPhysVol, &fHistory );
}
   */
// ********************************************************************
// LocateGlobalPointAndUpdateTouchable
// ********************************************************************
//
    /*
void GPNavigator_LocateGlobalPointAndUpdateTouchable(
                           const GPThreeVector&       position,
                                 G4VTouchable*        touchableToUpdate,
                           const G4bool               RelativeSearch )
{
  GPVPhysicalVolume* pPhysVol;
  pPhysVol = LocateGlobalPointAndSetup( position, 0, RelativeSearch);  
  touchableToUpdate->UpdateYourself( pPhysVol, &fHistory );
}
    */
// ********************************************************************
// GetVerboseLevel
// ********************************************************************
//

/*
G4int GPNavigator_GetVerboseLevel(GPNavigator *This)
{
  return This->fVerbose;
}
*/

// ********************************************************************
// SetVerboseLevel
// ********************************************************************
//

 /*
FQUALIFIER
void GPNavigator_SetVerboseLevel(GPNavigator *This, 
				 G4int level)
{
  This->fVerbose = level;
  GPNavigator_SetVerboseLevel(&(This->fnormalNav),level);
  //  fvoxelNav.SetVerboseLevel(level);
  //  fparamNav.SetVerboseLevel(level);
  //  freplicaNav.SetVerboseLevel(level);
  //  fregularNav.SetVerboseLevel(level);
}
 */

// ********************************************************************
// IsActive
// ********************************************************************
//
FQUALIFIER
G4bool GPNavigator_IsActive(GPNavigator *This)
{
  return This->fActive;
}

// ********************************************************************
// Activate
// ********************************************************************
//
FQUALIFIER
void GPNavigator_Activate(GPNavigator *This,
			  G4bool flag)
{
  This->fActive = flag;
}

// ********************************************************************
// EnteredDaughterVolume
//
// To inform the caller if the track is entering a daughter volume
// ********************************************************************
//
FQUALIFIER
G4bool GPNavigator_EnteredDaughterVolume(GPNavigator *This)
{
  return This->fEnteredDaughter;
}

// ********************************************************************
// ExitedMotherVolume
// ********************************************************************
//
FQUALIFIER
G4bool GPNavigator_ExitedMotherVolume(GPNavigator *This)
{
  return This->fExitedMother;
}

// ********************************************************************
// CheckMode
// ********************************************************************
//
/*
FQUALIFIER
void  GPNavigator_CheckMode(GPNavigator *This,
			    G4bool mode)
{
  This->fCheck = mode;
  GPNavigator_CheckMode(&(This->fnormalNav),mode);
//  fvoxelNav.CheckMode(mode);
//  fparamNav.CheckMode(mode);
//  freplicaNav.CheckMode(mode);
//  fregularNav.CheckMode(mode);
}
*/

// ********************************************************************
// IsCheckModeActive
// ********************************************************************
//
FQUALIFIER
G4bool GPNavigator_IsCheckModeActive(GPNavigator *This)
{
  return This->fCheck;
}

// ********************************************************************
// SetPushVerbosity
// ********************************************************************
//
FQUALIFIER
void GPNavigator_SetPushVerbosity(GPNavigator *This,
				  G4bool mode)
{
  This->fWarnPush = mode;
}

// ********************************************************************
// SeverityOfZeroStepping
//
// Reports on severity of error in case Navigator is stuck
// and is returning zero steps
// ********************************************************************
//
FQUALIFIER 
G4int GPNavigator_SeverityOfZeroStepping(GPNavigator *This,
					 G4int* noZeroSteps )
{
  G4int severity=0, noZeros= This->fNumberZeroSteps;
  if( noZeroSteps) *noZeroSteps = This->fNumberZeroSteps;

  if( noZeros >= This->fAbandonThreshold_NoZeroSteps )
  {
    severity = 10;
  }
  if( noZeros > 0 && noZeros < This->fActionThreshold_NoZeroSteps )
  {
    severity =  5 * noZeros / This->fActionThreshold_NoZeroSteps;
  }
  else if( noZeros == This->fActionThreshold_NoZeroSteps )
  {
    severity =  5; 
  }
  else if( noZeros >= This->fAbandonThreshold_NoZeroSteps - 2 )
  {
    severity =  9; 
  }
  else if( noZeros < This->fAbandonThreshold_NoZeroSteps - 2 )
  {
    severity =  5 + 4 * (noZeros- This->fAbandonThreshold_NoZeroSteps)
                      / This->fActionThreshold_NoZeroSteps;
  }
  return severity;
}

// ********************************************************************
// EnableBestSafety
// ********************************************************************
//
FQUALIFIER 
void GPNavigator_EnableBestSafety(GPNavigator *This,
				  G4bool value )
{
  //@@@G4FWO no voxel navigation
  ;
  //  fvoxelNav.EnableBestSafety( value );
}












