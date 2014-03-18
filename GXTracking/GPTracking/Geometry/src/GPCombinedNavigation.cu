
/* 
  These functions combine G4NormalNavigation and GPVoxelNavigation 
  G4FWP may not use this function, but keeps to retain Otto's original codes 
  - syjun
*/

#ifndef GPCOMBINEDNAVIGATION_CU
#define GPCOMBINEDNAVIGATION_CU

#include "GPCombinedNavigation.h"
#include "GPAuxiliaryNavServices.h"

// ********************************************************************
// ComputeStep
// ********************************************************************
//
FQUALIFIER 
G4double GPCombinedNavigation_ComputeStep(
		       GPVoxelNavigation *vox,
		       GPThreeVector localPoint,
		       GPThreeVector localDirection,
		       const G4double currentProposedStepLength,
		       G4double *newSafety,
		       GPNavigationHistory *history,
		       G4bool *validExitNormal,
		       GPThreeVector *exitNormal,
		       G4bool *exiting,
		       G4bool *entering,
		       GEOMETRYLOC GPVPhysicalVolume *(*pBlockedPhysical) )
{
  GEOMETRYLOC GPVPhysicalVolume *motherPhysical, *samplePhysical,
	*blockedExitedVol = GEOMETRYNULL;
  GEOMETRYLOC GPLogicalVolume *motherLogical;
  GEOMETRYLOC GPVSolid *motherSolid;
  GPThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int sampleNo;

  G4bool initialNode, noStep;
  GEOMETRYLOC GPSmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;

  motherPhysical = GPNavigationHistory_GetTopVolume( history );
  motherLogical = GPVPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = GPLogicalVolume_GetSolid(motherLogical);
  
  const G4bool voxelized = GPLogicalVolume_GetVoxelHeader(motherLogical) != GEOMETRYNULL;

  //
  // Compute mother safety
  //

  motherSafety = GPVSolid_DistanceToOut(motherSolid, localPoint);
  ourSafety = motherSafety;                 // Working isotropic safety

  //
  // Compute daughter safeties & intersections
  //

  // Exiting normal optimisation
  //
  if ( *exiting && *validExitNormal )
  {
    if ( GPThreeVector_dot(localDirection,*exitNormal)>=kMinExitingNormalCosine )
    {
      // Block exited daughter volume
      //
      blockedExitedVol = *pBlockedPhysical;
      ourSafety = 0;
    }
  }
  *exiting = false;
  *entering = false;

  G4int localNoDaughters = GPLogicalVolume_GetNoDaughters(motherLogical);
  
  initialNode = true;
  noStep = true;

  while (noStep)
  {
    if ( voxelized )
      {
	curVoxelNode = vox->fVoxelNode;
	curNoVolumes = GPSmartVoxelNode_GetNoContained(curVoxelNode);
      }
    else
      {
	curNoVolumes = localNoDaughters;
	noStep = false;
      }
    
    for (contentNo=curNoVolumes-1; contentNo>=0; contentNo--)
      {
	if (voxelized)
	  {
	    sampleNo = GPSmartVoxelNode_GetVolume(curVoxelNode,contentNo);
	  }
	else
	  {
	    sampleNo = contentNo;
	  }
	
	samplePhysical = GPLogicalVolume_GetDaughter(motherLogical,sampleNo);
        if ( samplePhysical!=blockedExitedVol )
	  {
	    GPAffineTransform sampleTf;
	    GPAffineTransform_Constructor3(&sampleTf,
					   GPVPhysicalVolume_GetRotation(samplePhysical),
					   GPVPhysicalVolume_GetTranslation(samplePhysical));
	    
	    GPAffineTransform_Invert(&sampleTf);
          
	    const GPThreeVector samplePoint =
	      GPAffineTransform_TransformPoint(&sampleTf,localPoint);
                     
	    GEOMETRYLOC const GPVSolid *sampleSolid =
	      GPLogicalVolume_GetSolid(GPVPhysicalVolume_GetLogicalVolume(samplePhysical ));
	    
	    const G4double sampleSafety =
	      GPVSolid_DistanceToIn(sampleSolid,samplePoint);

	    if ( sampleSafety<ourSafety )
	      {
		ourSafety = sampleSafety;
	      }
	    if ( sampleSafety<=ourStep )
	      {
		sampleDirection =
		  GPAffineTransform_TransformAxis( &sampleTf, localDirection );
		
		G4double sampleStep =
		  GPVSolid_DistanceToIn2(sampleSolid, samplePoint, sampleDirection);
		
		if ( sampleStep<=ourStep )
		  {
		    ourStep = sampleStep;
		    *entering = true;
		    *exiting = false;
		    *pBlockedPhysical = samplePhysical;
		    
		  }
	      }
	  }
      }
    if (initialNode)
      {
	if (voxelized)
	  {
	    initialNode = false;
	    voxelSafety = GPVoxelNavigation_ComputeVoxelSafety(vox,localPoint);
	    if ( voxelSafety<ourSafety )
	      {
		ourSafety = voxelSafety;
	      }
	  }
	if ( currentProposedStepLength<ourSafety )
	  {
	    // Guaranteed physics limited
	    //      
	    noStep = false;
	    *entering = false;
	    *exiting = false;
	    *pBlockedPhysical = GEOMETRYNULL;
	    ourStep = kInfinity;
	  }
	else
	  {
	    //
	    // Compute mother intersection if required
	    //
	    if ( motherSafety<=ourStep )
	      {
		G4double motherStep =
		  GPVSolid_DistanceToOut2( motherSolid, localPoint, localDirection,
					   true, validExitNormal, exitNormal);
		
		if ( motherStep<=ourStep )
		  {
		    ourStep = motherStep;
		    *exiting = true;
		    *entering = false;
		    if ( *validExitNormal )
		      {
			GPRotationMatrix rot = GPVPhysicalVolume_GetObjectRotationValue(motherPhysical);
			GPRotationMatrix inv = GPRotationMatrix_inverse(&rot);
			*exitNormal = GPRotationMatrix_apply( &inv, *exitNormal );
		      }
		  }
		else
		  {
		    *validExitNormal = false;
		  }
	      }
	  }
	*newSafety = ourSafety;
      }
    if ( noStep ) // not voxelized => noStep == false
      {
	noStep = GPVoxelNavigation_LocateNextVoxel(vox, localPoint, localDirection, ourStep);
      }
  }  // end -while (noStep)- loop
  
  return ourStep;
}


// ********************************************************************
// LevelLocate
// ********************************************************************
//
FQUALIFIER
G4bool GPCombinedNavigation_LevelLocate(
			GPVoxelNavigation *vox,
			GPNavigationHistory* history,
			GEOMETRYLOC const GPVPhysicalVolume* blockedVol,
			GPThreeVector globalPoint,
			GPThreeVector* globalDirection,
			const G4bool pLocatedOnEdge, 
			GPThreeVector *localPoint )
{
  GEOMETRYLOC GPSmartVoxelHeader *targetVoxelHeader;
  GEOMETRYLOC GPSmartVoxelNode *targetVoxelNode;
  GEOMETRYLOC GPVPhysicalVolume *targetPhysical, *samplePhysical;
  GEOMETRYLOC GPLogicalVolume *targetLogical;
  GEOMETRYLOC GPVSolid *sampleSolid;
  GPThreeVector samplePoint;
  G4int targetNoDaughters, sampleLogicalNo;
  
  targetPhysical = GPNavigationHistory_GetTopVolume(history);
  targetLogical = GPVPhysicalVolume_GetLogicalVolume(targetPhysical);
  targetVoxelHeader = GPLogicalVolume_GetVoxelHeader(targetLogical);
  
  const G4bool voxelized = targetVoxelHeader != GEOMETRYNULL;

  // Find the voxel containing the point
  //
  if ( voxelized )
  {
    targetVoxelNode =
      GPVoxelNavigation_VoxelLocate(vox,targetVoxelHeader,*localPoint);
    
    targetNoDaughters=GPSmartVoxelNode_GetNoContained(targetVoxelNode);
  }
  else {
    targetNoDaughters = GPLogicalVolume_GetNoDaughters(targetLogical);
  }
  
  if ( targetNoDaughters==0 ) return false;
  
  //
  // Search daughters in volume
  //
  for ( int sampleNo=targetNoDaughters-1; sampleNo>=0; sampleNo-- )
    {
      if ( voxelized )
	{
	  sampleLogicalNo = GPSmartVoxelNode_GetVolume(targetVoxelNode,sampleNo);
	}
      else
	{
	  sampleLogicalNo = sampleNo;
	}
	
      samplePhysical =
    	GPLogicalVolume_GetDaughter( targetLogical, sampleLogicalNo );
      
      if ( samplePhysical!=blockedVol )
	{
	  // Setup history
	  //
	  GPNavigationHistory_NewLevel(history, samplePhysical, kNormal, 0);
	  
	  sampleSolid =
	    GPLogicalVolume_GetSolid(
				     GPVPhysicalVolume_GetLogicalVolume( samplePhysical ));
	  
	  GPAffineTransform tf =
	    GPNavigationHistory_GetTopTransform( history );
	  
	  samplePoint =
	    GPAffineTransform_TransformPoint( &tf, globalPoint );
	  
	  if( GPAuxiliaryNavServices_CheckPointOnSurface(
							 sampleSolid, samplePoint, globalDirection, 
							 tf, pLocatedOnEdge) )
	    {
	      // Enter this daughter
	      //
	      *localPoint = samplePoint;
	      return true;
	    }
	  else
	    {
	      GPNavigationHistory_BackLevel( history );
	    }
	}
    }
  return false;
}

// ********************************************************************
// ComputeSafety
//
// Calculates the isotropic distance to the nearest boundary from the
// specified point in the local coordinate system. 
// The localpoint utilised must be within the current volume.
// ********************************************************************
//
FQUALIFIER 
G4double GPCombinedNavigation_ComputeSafety(
			GPVoxelNavigation *vox,
			GPThreeVector localPoint,
			GPNavigationHistory *history )
{	
  GEOMETRYLOC GPVPhysicalVolume *motherPhysical, *samplePhysical;
  GEOMETRYLOC GPLogicalVolume *motherLogical;
  GEOMETRYLOC GPVSolid *motherSolid;
  G4double motherSafety, ourSafety;
  G4int sampleNo;
  GEOMETRYLOC GPSmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;

  motherPhysical = GPNavigationHistory_GetTopVolume(history);
  motherLogical = GPVPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = GPLogicalVolume_GetSolid(motherLogical);
  const G4bool voxelized = GPLogicalVolume_GetVoxelHeader(motherLogical) != GEOMETRYNULL;

  //
  // Compute mother safety
  //

  motherSafety = GPVSolid_DistanceToOut(motherSolid,localPoint);
  ourSafety = motherSafety;                 // Working isotropic safety

  //
  // Compute daughter safeties 
  //

  if (voxelized)
  {
    //  Look only inside the current Voxel only (in the first version).
    //
    curVoxelNode = vox->fVoxelNode;
    curNoVolumes = GPSmartVoxelNode_GetNoContained(curVoxelNode);
  }
  else
  {
    curNoVolumes = GPLogicalVolume_GetNoDaughters(motherLogical);
  }

  for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
  {
    if ( voxelized )
      sampleNo = GPSmartVoxelNode_GetVolume(curVoxelNode,contentNo);
    else
      sampleNo = contentNo;
    
    samplePhysical = GPLogicalVolume_GetDaughter(motherLogical,sampleNo);
    
    GPAffineTransform sampleTf;
    GPAffineTransform_Constructor3(&sampleTf,
				   GPVPhysicalVolume_GetRotation(samplePhysical),
				   GPVPhysicalVolume_GetTranslation(samplePhysical));
	
    GPAffineTransform_Invert( &sampleTf );
    
    const GPThreeVector samplePoint =
      GPAffineTransform_TransformPoint(&sampleTf, localPoint);
    
    GEOMETRYLOC const GPVSolid *sampleSolid =
      GPLogicalVolume_GetSolid(
			       GPVPhysicalVolume_GetLogicalVolume( samplePhysical ));
    
    G4double sampleSafety =
      GPVSolid_DistanceToIn(sampleSolid, samplePoint);
    
    if ( sampleSafety<ourSafety )
      {
	ourSafety = sampleSafety;
      }
    
  }
  if ( voxelized )
    {
      voxelSafety = GPVoxelNavigation_ComputeVoxelSafety(vox,localPoint);
      if ( voxelSafety<ourSafety )
	{
	  ourSafety = voxelSafety;
	}
    }
  return ourSafety;
}

#endif
