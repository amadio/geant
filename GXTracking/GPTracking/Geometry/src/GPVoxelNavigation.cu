/**
 * G4VoxelNavigation implementation
 * based on G4VoxelNavigation.(i)cc of Geant 4.9.3
*/

#include "GPVoxelNavigation.h"
#include "GPNavigator.h"
#include "GPAuxiliaryNavServices.h"
#include "GPSmartVoxelNode.h"

// ********************************************************************
// Constructor
// ********************************************************************
//
FQUALIFIER
void GPVoxelNavigation_Constructor( GPVoxelNavigation *This )
{
  This->fVoxelDepth = -1;
  This->fVoxelNode = GEOMETRYNULL;
	
  //#ifdef USE_BLIST
  //  This->fBlist = NULL;
  //  This->fBlistSz = 0;
  //#endif
}

/*
#ifdef USE_BLIST
FQUALIFIER
void GPVoxelNavigation_EnlargeAndResetBlist( GPVoxelNavigation *This, G4int n )
{
  // TODO: knowingly leaking memory
  if ( This->fBlistSz < n ) {
    This->fBlist = realloc( This->fBlist, sizeof(char)*n );
    This->fBlistSz = n;
  }
  memset( This->fBlist, 0, sizeof(char)*n );
}
#endif
*/

#define std_min(a,b) (((a)<(b))?(a):(b))

// ********************************************************************
// ComputeVoxelSafety
//
// Computes safety from specified point to voxel boundaries
// using already located point
// o collected boundaries for most derived level
// o adjacent boundaries for previous levels
// ********************************************************************
//
FQUALIFIER
G4double GPVoxelNavigation_ComputeVoxelSafety(GPVoxelNavigation *This,
					      GPThreeVector localPoint)
{
  GEOMETRYLOC GPSmartVoxelHeader *curHeader;
  G4double voxelSafety, curNodeWidth;
  G4double curNodeOffset, minCurCommonDelta, maxCurCommonDelta;
  G4int minCurNodeNoDelta, maxCurNodeNoDelta;
  G4int localVoxelDepth, curNodeNo;
  EAxis curHeaderAxis;

  localVoxelDepth = This->fVoxelDepth;

  curHeader = This->fVoxelHeaderStack[localVoxelDepth];
  curHeaderAxis = This->fVoxelAxisStack[localVoxelDepth];
  curNodeNo = This->fVoxelNodeNoStack[localVoxelDepth];
  curNodeWidth = This->fVoxelSliceWidthStack[localVoxelDepth];
  
  // Compute linear intersection distance to boundaries of max/min
  // to collected nodes at current level
  //
  curNodeOffset = curNodeNo*curNodeWidth;
  
  maxCurNodeNoDelta = GPSmartVoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode)-curNodeNo;
  minCurNodeNoDelta = curNodeNo-GPSmartVoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode);
  minCurCommonDelta = GPThreeVector_coord(localPoint,curHeaderAxis)
	- GPSmartVoxelHeader_GetMinExtent(curHeader) - curNodeOffset;
	
  maxCurCommonDelta = curNodeWidth-minCurCommonDelta;

  if ( minCurNodeNoDelta<maxCurNodeNoDelta )
  {
    voxelSafety = minCurNodeNoDelta*curNodeWidth;
    voxelSafety += minCurCommonDelta;
  }
  else if (maxCurNodeNoDelta < minCurNodeNoDelta)
       {
         voxelSafety = maxCurNodeNoDelta*curNodeWidth;
         voxelSafety += maxCurCommonDelta;
        }
        else    // (maxCurNodeNoDelta == minCurNodeNoDelta)
        {
          voxelSafety = minCurNodeNoDelta*curNodeWidth;
          voxelSafety += std_min(minCurCommonDelta,maxCurCommonDelta);
        }

  // Compute isotropic safety to boundaries of previous levels
  // [NOT to collected boundaries]
  //
  while ( (localVoxelDepth>0) && (voxelSafety>0) )
  {
    localVoxelDepth--;
    curHeader = This->fVoxelHeaderStack[localVoxelDepth];
    curHeaderAxis = This->fVoxelAxisStack[localVoxelDepth];
    curNodeNo = This->fVoxelNodeNoStack[localVoxelDepth];
    curNodeWidth = This->fVoxelSliceWidthStack[localVoxelDepth];
    curNodeOffset = curNodeNo*curNodeWidth;
    minCurCommonDelta = GPThreeVector_coord(localPoint,curHeaderAxis)
                        - GPSmartVoxelHeader_GetMinExtent(curHeader) - curNodeOffset;
                        
    maxCurCommonDelta = curNodeWidth-minCurCommonDelta;
    
    if ( minCurCommonDelta<voxelSafety )
    {
      voxelSafety = minCurCommonDelta;
    }
    if ( maxCurCommonDelta<voxelSafety )
    {
      voxelSafety = maxCurCommonDelta;
    }
  }
  if ( voxelSafety<0 )
  {
    voxelSafety = 0;
  }

  return voxelSafety;
}

// ********************************************************************
// LevelLocate
// ********************************************************************
//
FQUALIFIER
G4bool GPVoxelNavigation_LevelLocate(GPVoxelNavigation *This,
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
  G4int targetNoDaughters;
  
  targetPhysical = GPNavigationHistory_GetTopVolume(history);
  targetLogical = GPVPhysicalVolume_GetLogicalVolume(targetPhysical);
  targetVoxelHeader = GPLogicalVolume_GetVoxelHeader(targetLogical);

  // Find the voxel containing the point
  //
  targetVoxelNode =
	GPVoxelNavigation_VoxelLocate(This,targetVoxelHeader,*localPoint);

  targetNoDaughters=GPSmartVoxelNode_GetNoContained(targetVoxelNode);
  if ( targetNoDaughters==0 ) return false;

  //
  // Search daughters in volume
  //
  for ( int sampleNo=targetNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    samplePhysical =
		GPLogicalVolume_GetDaughter( targetLogical, 
			GPSmartVoxelNode_GetVolume(targetVoxelNode,sampleNo));
                     
    if ( samplePhysical!=blockedVol )
    {
      // Setup history
      //
      GPNavigationHistory_NewLevel(history, samplePhysical, kNormal, 0);
      
      sampleSolid =
		GPLogicalVolume_GetSolid(
			GPVPhysicalVolume_GetLogicalVolume( samplePhysical ));
			
	  GPAffineTransform tf = GPNavigationHistory_GetTopTransform( history );
			
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
  
  //return G4NormalNavigation_LevelLocate( history, blockedVol, &globalPoint, globalDirection, pLocatedOnEdge, localPoint );
}

// ********************************************************************
// ComputeStep
// ********************************************************************
//
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
				       GEOMETRYLOC GPVPhysicalVolume *(*pBlockedPhysical) )
{
  /*return G4NormalNavigation_ComputeStep(
		localPoint, localDirection, currentProposedStepLength,
		newSafety, history, validExitNormal, exitNormal,
		exiting, entering, pBlockedPhysical );*/
		
  GEOMETRYLOC GPVPhysicalVolume *motherPhysical, *samplePhysical, *blockedExitedVol = GEOMETRYNULL ;
  GEOMETRYLOC GPLogicalVolume *motherLogical;
  GEOMETRYLOC GPVSolid *motherSolid;
  GPThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int sampleNo; // , localNoDaughters;

  G4bool initialNode, noStep;
  GEOMETRYLOC GPSmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;

  motherPhysical = GPNavigationHistory_GetTopVolume( history );
  motherLogical = GPVPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = GPLogicalVolume_GetSolid(motherLogical);

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

  //#ifdef USE_BLIST
  //  GPVoxelNavigation_EnlargeAndResetBlist( This, GPLogicalVolume_GetNoDaughters(motherLogical) );
  //#endif

  initialNode = true;
  noStep = true;

  while (noStep)
  {
    curVoxelNode = This->fVoxelNode;
    curNoVolumes = GPSmartVoxelNode_GetNoContained(curVoxelNode);
    for (contentNo=curNoVolumes-1; contentNo>=0; contentNo--)
    {
      sampleNo = GPSmartVoxelNode_GetVolume(curVoxelNode,contentNo);
      
      //#ifdef USE_BLIST
      //      if (!This->fBlist[sampleNo])
      //      {
      //		This->fBlist[sampleNo] = 1;
      //#endif
		
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
				GPLogicalVolume_GetSolid(
					GPVPhysicalVolume_GetLogicalVolume(
						samplePhysical ));
          
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
	    //#ifdef USE_BLIST
	    //          } // -- FBLIST
	    //#endif
        }
      }
    }
    if (initialNode)
    {	
      initialNode = false;
      voxelSafety = GPVoxelNavigation_ComputeVoxelSafety(This,localPoint);
      if ( voxelSafety<ourSafety )
      {
        ourSafety = voxelSafety;
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
          G4double motherStep = GPVSolid_DistanceToOut2( motherSolid, 
							     localPoint, 
							     localDirection,
							     true, 
							     validExitNormal, 
							     exitNormal);

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
    if (noStep)
    {
      noStep = GPVoxelNavigation_LocateNextVoxel(This, localPoint, 
						 localDirection, ourStep);
    }
  }  // end -while (noStep)- loop

  return ourStep;
}

// class G4VoxelNavigation Inline implementation
//
// --------------------------------------------------------------------

// ********************************************************************
// VoxelLocate
// ********************************************************************
//
FQUALIFIER
GEOMETRYLOC GPSmartVoxelNode* GPVoxelNavigation_VoxelLocate(GPVoxelNavigation *This,
						GEOMETRYLOC GPSmartVoxelHeader* pHead,
						GPThreeVector localPoint)
{
  GEOMETRYLOC GPSmartVoxelHeader *targetVoxelHeader=pHead;
  GEOMETRYLOC GPSmartVoxelNode *targetVoxelNode = GEOMETRYNULL;
  GEOMETRYLOC GPSmartVoxelProxy *sampleProxy;
  EAxis targetHeaderAxis;
  G4double targetHeaderMin, targetHeaderNodeWidth;
  G4int targetHeaderNoSlices, targetNodeNo;

  This->fVoxelDepth = 0;

  //  while ( !targetVoxelNode )
  while ( targetVoxelNode == GEOMETRYNULL )
  {
    targetHeaderAxis = GPSmartVoxelHeader_GetAxis(targetVoxelHeader);
    targetHeaderNoSlices = GPSmartVoxelHeader_GetNoSlices(targetVoxelHeader);
    targetHeaderMin = GPSmartVoxelHeader_GetMinExtent(targetVoxelHeader);
    targetHeaderNodeWidth =
		(GPSmartVoxelHeader_GetMaxExtent(targetVoxelHeader)-targetHeaderMin)
                          / targetHeaderNoSlices;
    targetNodeNo = (G4int)(
      (GPThreeVector_coord(localPoint,targetHeaderAxis)-targetHeaderMin)
                          / targetHeaderNodeWidth);
                          
    // Rounding protection
    //
    if ( targetNodeNo<0 )
    {
		targetNodeNo = 0;
    }
    else if ( targetNodeNo>=targetHeaderNoSlices )
	{
		targetNodeNo = targetHeaderNoSlices-1;
	}
         
    // Stack info for stepping
    //
    
    This->fVoxelAxisStack[This->fVoxelDepth] = targetHeaderAxis;
    This->fVoxelNoSlicesStack[This->fVoxelDepth] = targetHeaderNoSlices;
    This->fVoxelSliceWidthStack[This->fVoxelDepth] = targetHeaderNodeWidth;
    This->fVoxelNodeNoStack[This->fVoxelDepth] = targetNodeNo;
    This->fVoxelHeaderStack[This->fVoxelDepth] = targetVoxelHeader;
    sampleProxy = GPSmartVoxelHeader_GetSlice(targetVoxelHeader, targetNodeNo);

    if ( GPSmartVoxelProxy_IsNode(sampleProxy) )
    {
      targetVoxelNode = GPSmartVoxelProxy_GetNode(sampleProxy);
    }
    else
    {
      targetVoxelHeader = GPSmartVoxelProxy_GetHeader(sampleProxy);
      This->fVoxelDepth++;
      //      myAssert(This->fVoxelDepth < K_MAX_VOXEL_STACK_DEPTH);
    }
  }
  
  This->fVoxelNode = targetVoxelNode;
  return targetVoxelNode;
}


// ********************************************************************
// LocateNextVoxel
//
// Finds the next voxel from the current voxel and point
// in the specified direction
//
// Returns false if all voxels considered
//              [current Step ends inside same voxel or leaves all voxels]
//         true  otherwise
//              [the information on the next voxel is put into the set of
//               fVoxel* variables & "stacks"] 
// ********************************************************************
// 
FQUALIFIER
G4bool GPVoxelNavigation_LocateNextVoxel(GPVoxelNavigation *This,
					 GPThreeVector localPoint,
					 GPThreeVector localDirection,
					 const G4double currentStep )
{
  GEOMETRYLOC GPSmartVoxelHeader *workHeader= GEOMETRYNULL, *newHeader= GEOMETRYNULL;
  GEOMETRYLOC GPSmartVoxelProxy *newProxy= GEOMETRYNULL;
  GEOMETRYLOC GPSmartVoxelNode *newVoxelNode= GEOMETRYNULL;
  GPThreeVector targetPoint, voxelPoint;
  G4double workNodeWidth, workMinExtent, workCoord;
  G4double minVal, maxVal, newDistance=0.;
  G4double newHeaderMin, newHeaderNodeWidth;
  G4int depth=0, newDepth=0, workNodeNo=0, newNodeNo=0, newHeaderNoSlices=0;
  EAxis workHeaderAxis, newHeaderAxis;
  G4bool isNewVoxel=false;
  
  G4double currentDistance = currentStep;

  // Determine if end of Step within current voxel
  //
  for (depth=0; depth<This->fVoxelDepth; depth++)
  {
    targetPoint =
		GPThreeVector_saxpy(currentDistance,localDirection,localPoint);
			
    newDistance = currentDistance;
    workHeader = This->fVoxelHeaderStack[depth];
    workHeaderAxis = This->fVoxelAxisStack[depth];
    workNodeNo = This->fVoxelNodeNoStack[depth];
    workNodeWidth = This->fVoxelSliceWidthStack[depth];
    workMinExtent = GPSmartVoxelHeader_GetMinExtent(workHeader);
    workCoord = GPThreeVector_coord(targetPoint,workHeaderAxis);
    minVal = workMinExtent+workNodeNo*workNodeWidth;

    if ( minVal<=workCoord+kCarTolerance*0.5 )
    {
      maxVal = minVal+workNodeWidth;
      if ( maxVal<=workCoord-kCarTolerance*0.5 )
      {
        // Must consider next voxel
        //
        newNodeNo = workNodeNo+1;
        newHeader = workHeader;
        newDistance = (maxVal-GPThreeVector_coord(localPoint,workHeaderAxis))
                    / GPThreeVector_coord(localDirection,workHeaderAxis);
        isNewVoxel = true;
        newDepth = depth;
      }
    }
    else
    {
      newNodeNo = workNodeNo-1;
      newHeader = workHeader;
      newDistance = (minVal-GPThreeVector_coord(localPoint,workHeaderAxis))
                  / GPThreeVector_coord(localDirection,workHeaderAxis);
      isNewVoxel = true;
      newDepth = depth;
    }
    currentDistance = newDistance;
  }
  targetPoint = 
	GPThreeVector_saxpy(currentDistance,localDirection,localPoint);

  // Check if end of Step within collected boundaries of current voxel
  //
  depth = This->fVoxelDepth;
  {
    workHeader = This->fVoxelHeaderStack[depth];
    workHeaderAxis = This->fVoxelAxisStack[depth];
    workNodeNo = This->fVoxelNodeNoStack[depth];
    workNodeWidth = This->fVoxelSliceWidthStack[depth];
    workMinExtent = GPSmartVoxelHeader_GetMinExtent(workHeader);
    workCoord = GPThreeVector_coord(targetPoint,workHeaderAxis);
    minVal = workMinExtent+GPSmartVoxelNode_GetMinEquivalentSliceNo(This->fVoxelNode)*workNodeWidth;

    if ( minVal<=workCoord+kCarTolerance*0.5 )
    {
      maxVal = workMinExtent+(GPSmartVoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode)+1)
                            *workNodeWidth;
      if ( maxVal<=workCoord-kCarTolerance*0.5 )
      {
        newNodeNo = GPSmartVoxelNode_GetMaxEquivalentSliceNo(This->fVoxelNode)+1;
        newHeader = workHeader;
        newDistance = (maxVal-GPThreeVector_coord(localPoint,workHeaderAxis))
                    / GPThreeVector_coord(localDirection,workHeaderAxis);
        isNewVoxel = true;
        newDepth = depth;
      }
    }
    else
    {
      newNodeNo = GPSmartVoxelNode_GetMinEquivalentSliceNo(This->fVoxelNode)-1;
      newHeader = workHeader;
      newDistance = (minVal-GPThreeVector_coord(localPoint,workHeaderAxis))
                  / GPThreeVector_coord(localDirection,workHeaderAxis);
      isNewVoxel = true;
      newDepth = depth;
    }
    currentDistance = newDistance;
  }
  if (isNewVoxel)
  {
    // Compute new voxel & adjust voxel stack
    //
    // newNodeNo=Candidate node no at 
    // newDepth =refinement depth of crossed voxel boundary
    // newHeader=Header for crossed voxel
    // newDistance=distance to crossed voxel boundary (along the track)
    //
    if ( (newNodeNo<0) || (newNodeNo>=GPSmartVoxelHeader_GetNoSlices(newHeader)))
    {
      // Leaving mother volume
      //
      isNewVoxel = false;
    }
    else
    {
      // Compute intersection point on the least refined
      // voxel boundary that is hit
      //
      voxelPoint = GPThreeVector_saxpy(newDistance,localDirection,localPoint);
      //      myAssert(newDepth < K_MAX_VOXEL_STACK_DEPTH);
      This->fVoxelNodeNoStack[newDepth] = newNodeNo;
      This->fVoxelDepth = newDepth;
      newVoxelNode = 0;
      //      while ( !newVoxelNode )
      while ( newVoxelNode == GEOMETRYNULL )
      {
        newProxy = GPSmartVoxelHeader_GetSlice(newHeader,newNodeNo);
        if ( GPSmartVoxelProxy_IsNode(newProxy) )
        {
          newVoxelNode = GPSmartVoxelProxy_GetNode(newProxy);
        }
        else
        {
          This->fVoxelDepth++;
	  //          myAssert(This->fVoxelDepth < K_MAX_VOXEL_STACK_DEPTH);
          newHeader = GPSmartVoxelProxy_GetHeader(newProxy);
          newHeaderAxis = GPSmartVoxelHeader_GetAxis(newHeader);
          newHeaderNoSlices = GPSmartVoxelHeader_GetNoSlices(newHeader);
          newHeaderMin = GPSmartVoxelHeader_GetMinExtent(newHeader);
          newHeaderNodeWidth =
			(GPSmartVoxelHeader_GetMaxExtent(newHeader)-newHeaderMin)
                             / newHeaderNoSlices;
          newNodeNo = (G4int)(
			(GPThreeVector_coord(voxelPoint,newHeaderAxis)-newHeaderMin)
                             / newHeaderNodeWidth );
          // Rounding protection
          //
          if ( newNodeNo<0 )
          {
            newNodeNo=0;
          }
          else if ( newNodeNo>=newHeaderNoSlices )
               {
                 newNodeNo = newHeaderNoSlices-1;
               }
          // Stack info for stepping
          //
          This->fVoxelAxisStack[This->fVoxelDepth] = newHeaderAxis;
          This->fVoxelNoSlicesStack[This->fVoxelDepth] = newHeaderNoSlices;
          This->fVoxelSliceWidthStack[This->fVoxelDepth] = newHeaderNodeWidth;
          This->fVoxelNodeNoStack[This->fVoxelDepth] = newNodeNo;
          This->fVoxelHeaderStack[This->fVoxelDepth] = newHeader;
        }
      }
      This->fVoxelNode = newVoxelNode;
    }
  }
  return isNewVoxel;        
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
G4double GPVoxelNavigation_ComputeSafety(GPVoxelNavigation *This,
					 GPThreeVector localPoint,
					 GPNavigationHistory *history ,
					 const G4double       maxLength)
{
  GEOMETRYLOC GPVPhysicalVolume *motherPhysical, *samplePhysical;
  GEOMETRYLOC GPLogicalVolume *motherLogical;
  GEOMETRYLOC GPVSolid *motherSolid;
  G4double motherSafety, ourSafety;
  G4int sampleNo; // , localNoDaughters;
  GEOMETRYLOC GPSmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
  G4double voxelSafety;

  motherPhysical = GPNavigationHistory_GetTopVolume(history);
  motherLogical = GPVPhysicalVolume_GetLogicalVolume(motherPhysical);
  motherSolid = GPLogicalVolume_GetSolid(motherLogical);

  //
  // Compute mother safety
  //

  motherSafety = GPVSolid_DistanceToOut(motherSolid,localPoint);
  ourSafety = motherSafety;                 // Working isotropic safety

  //
  // Compute daughter safeties 
  //

  //localNoDaughters = GPLogicalVolume_GetNoDaughters(motherLogical);

  //  Look only inside the current Voxel only (in the first version).
  //
  curVoxelNode = This->fVoxelNode;
  curNoVolumes = GPSmartVoxelNode_GetNoContained(curVoxelNode);

  for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
  {
    sampleNo = GPSmartVoxelNode_GetVolume(curVoxelNode,contentNo);
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
  voxelSafety = GPVoxelNavigation_ComputeVoxelSafety(This,localPoint);
  if ( voxelSafety<ourSafety )
  {
    ourSafety = voxelSafety;
  }
  return ourSafety;
}
