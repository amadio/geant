/**
 * G4VoxelHeader (host side) implementation
 * based on G4VoxelHeader.cc of Geant 4.9.3
 */

#ifndef GPSMARTVOXELHEADER_CPP
#define GPSMARTVOXELHEADER_CPP

//#include "G4BuildVoxels.h"
/*
#include "GPVoxeldefs.h"
#include "GPVoxelLimits.h"
#include "GPLogicalVolume.h"
#include "GPVPhysicalVolume.h"
#include "GPAffineTransform.h"

#include "GPSmartVoxelHeader.h"
#include "GPSmartVoxelNode.h"
#include "GPSmartVoxelProxy.h"
*/
//#include "GPLogicalVolume_inline.c"
//#include "GPVSolid_inline.c"

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <assert.h>
#include <cstring>

#include "GPVoxelHeader.h"

void GPSmartVoxelHeader_SetMinVoxelLimits( int lev1, int lev2, int lev3 )
{
  minVoxelVolumesLevel1 = lev1;
  minVoxelVolumesLevel2 = lev2;
  minVoxelVolumesLevel3 = lev3;
}

int GPSmartVoxelHeader_GetMinVoxelVolumesLevel1()
{
  return minVoxelVolumesLevel1;
}

int GPSmartVoxelHeader_GetMinVoxelVolumesLevel2()
{
  return minVoxelVolumesLevel2;
}

int GPSmartVoxelHeader_GetMinVoxelVolumesLevel3()
{
  return minVoxelVolumesLevel3;
}

// ***************************************************************************
// Constructor for topmost header, to begin voxel construction at a
// given logical volume.
// Constructs target List of volumes, calls "Build and refine" constructor.
// Assumes all daughters represent single volumes (ie. no divisions
// or parametric)
// ***************************************************************************
//
void GPSmartVoxelHeader_Constructor(GPSmartVoxelHeader *This,
				    GPLogicalVolume* pVolume,
				    G4int pSlice,
				    G4double smartless )
{
  assert( pVolume != NULL );
  
  This->fnumSlices = 0;
  This->fslices = NULL;
  
  This->faxis = kUndefined;
  This->fminEquivalent = pSlice;
  This->fmaxEquivalent = pSlice;
  This->fparamAxis = kUndefined;
  
  if ( GPLogicalVolume_GetNoDaughters(pVolume) >= minVoxelVolumesLevel1 ) {
    GPSmartVoxelHeader_BuildVoxels(This,pVolume,smartless);
  }
}

// ***************************************************************************
// Prs and refines voxels between specified limits, considering only
// the physical volumes numbered `pCandidates'. `pSlice' is used to set max
// and min equivalent slice nos for the header - they apply to the level
// of the header, not its nodes.
// ***************************************************************************
//
void GPSmartVoxelHeader_Constructor2(GPSmartVoxelHeader *This,
				     GPLogicalVolume* pVolume,
				     GPVoxelLimits pLimits,
				     const GPVolumeNosVector* pCandidates,
				     G4int pSlice,
				     G4double smartless)
{
  This->fnumSlices = 0;
  This->fslices = NULL;
  
  This->faxis = kUndefined;
  This->fminEquivalent = pSlice;
  This->fmaxEquivalent = pSlice;
  This->fparamAxis = kUndefined;
  
  GPSmartVoxelHeader_BuildVoxelsWithinLimits(This,pVolume,pLimits,pCandidates, smartless);
}

// ***************************************************************************
// Destructor:
// deletes all proxies and underlying objects.
// ***************************************************************************
//
void
GPSmartVoxelHeader_Destructor(GPSmartVoxelHeader *This)
{
  // Manually destroy underlying nodes/headers
  // Delete collected headers and nodes once only
  //
  G4int node, proxy, maxNode=GPSmartVoxelHeader_GetNoSlices(This);
  GPSmartVoxelProxy *lastProxy=0;
  GPSmartVoxelNode *dyingNode, *lastNode=0;
  GPSmartVoxelHeader *dyingHeader, *lastHeader=0;

  for (node=0; node<maxNode; node++)
  {
    if (GPSmartVoxelProxy_IsHeader(This->fslices[node]))
    {
      dyingHeader = GPSmartVoxelProxy_GetHeader(This->fslices[node]);
      if (lastHeader!=dyingHeader)
      {
        lastHeader = dyingHeader;
        lastNode = 0;
        GPSmartVoxelHeader_Destructor( dyingHeader );
        delete dyingHeader;
      }
    }
    else
    {
      dyingNode = GPSmartVoxelProxy_GetNode(This->fslices[node]);
      if (dyingNode!=lastNode)
      {
        lastNode=dyingNode;
        lastHeader=0;
        
        GPSmartVoxelNode_Destructor( dyingNode );
        delete dyingNode;
      }
    }
  }
  // Delete proxies
  //
  for (proxy=0; proxy<maxNode; proxy++)
  {
    if (This->fslices[proxy]!=lastProxy)
    {
      lastProxy = This->fslices[proxy];
      delete lastProxy;
    }
  }
  
  std::free(This->fslices);
}

// ***************************************************************************
// Equality operator: returns true if contents are equivalent.
// Implies a deep search through contained nodes/header.
// Compares headers' axes,sizes,extents. Returns false if different.
// For each contained proxy, determines whether node/header, compares and
// returns if different. Compares and returns if proxied nodes/headers
// are different.
// ***************************************************************************
//
G4bool GPSmartVoxelHeader_operator_equal(GPSmartVoxelHeader *This,
					 GPSmartVoxelHeader* pHead)
{
  if ( (GPSmartVoxelHeader_GetAxis(This)      == GPSmartVoxelHeader_GetAxis(pHead))
    && (GPSmartVoxelHeader_GetNoSlices(This)  == GPSmartVoxelHeader_GetNoSlices(This))
    && (GPSmartVoxelHeader_GetMinExtent(This) == GPSmartVoxelHeader_GetMinExtent(This))
    && (GPSmartVoxelHeader_GetMaxExtent(This) == GPSmartVoxelHeader_GetMaxExtent(This)) )
  {
    G4int node, maxNode;
    GPSmartVoxelProxy *leftProxy, *rightProxy;
    GPSmartVoxelHeader *leftHeader, *rightHeader;
    GPSmartVoxelNode *leftNode, *rightNode;

    maxNode=GPSmartVoxelHeader_GetNoSlices(This);
    for (node=0; node<maxNode; node++)
    {
      leftProxy  = GPSmartVoxelHeader_GetSlice(This,node);
      rightProxy = GPSmartVoxelHeader_GetSlice(pHead,node);
      if (GPSmartVoxelProxy_IsHeader(leftProxy))
      {
        if (GPSmartVoxelProxy_IsNode(rightProxy))
        {
          return false;
        }
        else
        {
          leftHeader  = GPSmartVoxelProxy_GetHeader(leftProxy);
          rightHeader = GPSmartVoxelProxy_GetHeader(rightProxy);
          if (!(GPSmartVoxelHeader_operator_equal(leftHeader,rightHeader)))
          {
            return false;
          }
        }
      }
      else
      {
        if (GPSmartVoxelProxy_IsHeader(rightProxy))
        {
          return false;
        }
        else
        {
          leftNode  = GPSmartVoxelProxy_GetNode(leftProxy);
          rightNode = GPSmartVoxelProxy_GetNode(rightProxy);
          if (!(GPSmartVoxelNode_operator_equal(leftNode,rightNode)))
          {
            return false;
          }
        }
      }
    }
    return true;
  }
  else
  {
    return false;
  }
}

void GPSmartVoxelHeader_SetSlices(GPSmartVoxelHeader *This,
				  const GPProxyVector *slices )
{
  This->fnumSlices = slices->size();
	
  std::size_t len = sizeof(GPSmartVoxelProxy*)*This->fnumSlices;
  
  This->fslices = (GPSmartVoxelProxy**)std::realloc( This->fslices, len );
  memcpy( This->fslices, &((*slices)[0]), len );
}

// ***************************************************************************
// Builds the nodes corresponding to slices between the specified limits
// and along the specified axis, using candidate volume no.s in the vector
// pCandidates. If the `daughters' are replicated volumes (ie. the logical
// volume has a single replicated/parameterised volume for a daughter)
// the candidate no.s are interpreted as PARAMETERISED volume no.s & 
// PARAMETERISATIONs are applied to compute transformations & solid
// dimensions appropriately. The volume must be parameterised - ie. has a
// parameterisation object & non-consuming) - in this case.
// 
// Returns pointer to built node "structure" (guaranteed non NULL) consisting
// of GPSmartVoxelNodeProxies referring to GPSmartVoxelNodes.
// ***************************************************************************
//
GPProxyVector* GPSmartVoxelHeader_BuildNodes(GPSmartVoxelHeader *This,
					     GPLogicalVolume* pVolume,
					     GPVoxelLimits pLimits,
					     const GPVolumeNosVector* pCandidates,
					     EAxis pAxis,
					     G4double smartlessUser )
{
  (void)This;
	
  G4double motherMinExtent= kInfinity, motherMaxExtent= -kInfinity,
           targetMinExtent= kInfinity, targetMaxExtent= -kInfinity;
  GPVPhysicalVolume *pDaughter=0;
  GPVSolid *targetSolid;
  GPAffineTransform targetTransform = GPAffineTransform_create_id();
  G4int nCandidates = pCandidates->size();
  G4int nVol, nNode, targetVolNo;
  GPVoxelLimits noLimits; GPVoxelLimits_Constructor(&noLimits);
    
  // Compute extent of logical volume's solid along this axis
  // NOTE: results stored locally and not preserved/reused
  //
  GPVSolid* outerSolid = GPLogicalVolume_GetSolid(pVolume);
  assert( outerSolid != NULL );
  
  GPAffineTransform origin = GPAffineTransform_create_id();
  if( !GPVSolid_CalculateExtent(outerSolid, pAxis, pLimits, origin,
                                   &motherMinExtent, &motherMaxExtent) )
  {
    GPVSolid_CalculateExtent(outerSolid, pAxis, noLimits, origin,
                                &motherMinExtent, &motherMaxExtent);
  }
  GPVolumeExtentVector minExtents(nCandidates,0.);
  GPVolumeExtentVector maxExtents(nCandidates,0.);
    
  // Compute extents
  //
  for (nVol=0; nVol<nCandidates; nVol++) {
    targetVolNo=(*pCandidates)[nVol];
    
    pDaughter=GPLogicalVolume_GetDaughter(pVolume,targetVolNo);
    
    // Setup daughter's transformations
    //
    GPAffineTransform targetTransform;
    GPAffineTransform_Constructor3(&targetTransform,
				   GPVPhysicalVolume_GetRotation(pDaughter),
				   GPVPhysicalVolume_GetTranslation(pDaughter));
    
    // Get underlying (and setup) solid
    //
    targetSolid = GPLogicalVolume_GetSolid(
		  GPVPhysicalVolume_GetLogicalVolume(pDaughter) );
    
    assert( targetSolid != NULL );
    
    // Calculate extents
    //
    if(!GPVSolid_CalculateExtent(targetSolid, pAxis, pLimits, targetTransform,
				 &targetMinExtent, &targetMaxExtent))
      {
	GPVSolid_CalculateExtent(targetSolid, pAxis, noLimits, targetTransform,
				 &targetMinExtent, &targetMaxExtent);
      }
    minExtents[nVol] = targetMinExtent;
    maxExtents[nVol] = targetMaxExtent;
    
    // Check not entirely outside mother when processing toplevel nodes
    //
    if( (!GPVoxelLimits_IsLimited(&pLimits)) && ((targetMaxExtent<=motherMinExtent) ||(targetMinExtent>=motherMaxExtent)) )
      {
	return NULL;
      }
    
  }

  // Extents of all daughters known

  // Calculate minimum slice width, only including volumes inside the limits
  //
  G4double minWidth = kInfinity;
  G4double currentWidth;
  for (nVol=0; nVol<nCandidates; nVol++)
  {
    // currentWidth should -always- be a positive value. Inaccurate computed extent
    // from the solid or situations of malformed geometries (overlaps) may lead to
    // negative values and therefore unpredictable crashes !
    //
    currentWidth = abs(maxExtents[nVol]-minExtents[nVol]);
    if ( (currentWidth<minWidth)
      && (maxExtents[nVol]>=GPVoxelLimits_GetMinExtent(&pLimits,pAxis))
      && (minExtents[nVol]<=GPVoxelLimits_GetMaxExtent(&pLimits,pAxis)) )
    {
      minWidth = currentWidth;
    }
  }

  // No. of Nodes formula - nearest integer to
  // mother width/half min daughter width +1
  //
  G4double noNodesExactD = ((motherMaxExtent-motherMinExtent)*2.0/minWidth)+1.0;

  // Compare with "smartless quality", i.e. the average number of slices
  // used per contained volume.
  //
  G4double smartlessComputed = noNodesExactD / nCandidates;
  
  G4double smartless = (smartlessComputed <= smartlessUser)
                       ? smartlessComputed : smartlessUser;
  
  G4double noNodesSmart = smartless*nCandidates;
  G4int    noNodesExactI = G4int(noNodesSmart);
  G4int    noNodes = ((noNodesSmart-noNodesExactI)>=0.5)
                     ? noNodesExactI+1 : noNodesExactI;
  if( noNodes == 0 ) { noNodes=1; }
  
  //  const G4int kMaxVoxelNodes = K_MAX_VOXEL_NODES;

  if (noNodes > kMaxVoxelNodes)
  {
    noNodes=kMaxVoxelNodes; // TODO?
  }
  G4double nodeWidth = (motherMaxExtent-motherMinExtent)/noNodes;

// Create GPSmartVoxelNodes. Will Add proxies before setting fslices
//
  GPNodeVector* nodeList = new GPNodeVector();
  nodeList->reserve(noNodes);
  
  for (nNode=0; nNode<noNodes; nNode++)
  {
    GPSmartVoxelNode *pNode;
    pNode = new GPSmartVoxelNode;
    GPSmartVoxelNode_Constructor(pNode,nNode);
    
    nodeList->push_back(pNode);
  }

  // All nodes created (empty)

  // Fill nodes: Step through extent lists
  //
  for (nVol=0; nVol<nCandidates; nVol++)
  {
    G4int nodeNo, minContainingNode, maxContainingNode;
    minContainingNode = G4int((minExtents[nVol]-motherMinExtent)/nodeWidth);
    maxContainingNode = G4int((maxExtents[nVol]-motherMinExtent)/nodeWidth);

    // Only add nodes that are inside the limits of the axis
    //
    if ( (maxContainingNode>=0) && (minContainingNode<noNodes) )
    {
      // If max extent is on max boundary => maxContainingNode=noNodes:
      // should be one less as nodeList has noNodes entries
      //
      if (maxContainingNode>=noNodes)
      {
        maxContainingNode = noNodes-1;
      }
      //
      // Protection against protruding volumes
      //
      if (minContainingNode<0)
      {
        minContainingNode=0;
      }
      for (nodeNo=minContainingNode; nodeNo<=maxContainingNode; nodeNo++)
      {
		  GPSmartVoxelNode_Insert( (*nodeList)[nodeNo], (*pCandidates)[nVol] );
      }
    }
  }

  // All nodes filled

  // Create proxy List : caller has deletion responsibility
  // (but we must delete nodeList *itself* - not the contents)
  //
  GPProxyVector* proxyList = new GPProxyVector();
  proxyList->reserve(noNodes);
  
  //
  // Fill proxy List
  //
  for (nNode=0; nNode<noNodes; nNode++)
  {

    GPSmartVoxelProxy* pProxyNode = new GPSmartVoxelProxy;
    GPSmartVoxelProxy_SetNode( pProxyNode, ((*nodeList)[nNode]) );
    
    proxyList->push_back(pProxyNode);
  }
  delete nodeList;
  return proxyList;
}


// ***************************************************************************
// Collects common nodes at our level, deleting all but one to save
// memory, and adjusting stored slice pointers appropriately.
//
// Preconditions:
// o the slices have not previously be "collected"
// o all of the slices are nodes.
// ***************************************************************************
//
void GPSmartVoxelHeader_CollectEquivalentNodes( GPSmartVoxelHeader *This )
{
  G4int sliceNo, maxNo, equivNo;
  G4int maxNode=GPSmartVoxelHeader_GetNoSlices(This);
  GPSmartVoxelNode *equivNode;
  GPSmartVoxelProxy *equivProxy;

  for (sliceNo=0; sliceNo<maxNode; sliceNo++)
  {
    equivProxy=This->fslices[sliceNo];

    // Assumption (see preconditions): all slices are nodes
    //
    equivNode = GPSmartVoxelProxy_GetNode(equivProxy);
    maxNo = GPSmartVoxelNode_GetMaxEquivalentSliceNo(equivNode);
    if (maxNo != sliceNo)
    {
      // Do collection between sliceNo and maxNo inclusive
      //
      for (equivNo=sliceNo+1; equivNo<=maxNo; equivNo++)
      {
		GPSmartVoxelNode *dyingNode = GPSmartVoxelProxy_GetNode(This->fslices[equivNo]);
		GPSmartVoxelNode_Destructor( dyingNode );
        delete dyingNode;
        delete This->fslices[equivNo];
        This->fslices[equivNo] = equivProxy;
      }
      sliceNo = maxNo;
    }
  }
}

// ***************************************************************************
// Collects common headers at our level, deleting all but one to save
// memory, and adjusting stored slice pointers appropriately.
// 
// Preconditions:
// o if a header forms part of a range of equivalent slices
//   (ie. GetMaxEquivalentSliceNo()>GetMinEquivalentSliceNo()),
//   it is assumed that all slices in the range are headers.
// o this will be true if a constant Expression is used to evaluate
//   when to refine nodes.
// ***************************************************************************
//
void GPSmartVoxelHeader_CollectEquivalentHeaders( GPSmartVoxelHeader *This )
{
  G4int sliceNo, maxNo, equivNo;
  G4int maxNode = GPSmartVoxelHeader_GetNoSlices(This);
  GPSmartVoxelHeader *equivHeader, *sampleHeader;
  GPSmartVoxelProxy *equivProxy;

  for (sliceNo=0; sliceNo<maxNode; sliceNo++)
  {
    equivProxy = This->fslices[sliceNo];
    if (GPSmartVoxelProxy_IsHeader(equivProxy))
    {
      equivHeader = GPSmartVoxelProxy_GetHeader(equivProxy);
      maxNo = GPSmartVoxelHeader_GetMaxEquivalentSliceNo(equivHeader);
      if (maxNo != sliceNo)
      {
        // Attempt collection between sliceNo and maxNo inclusive:
        // look for common headers. All slices between sliceNo and maxNo
        // are guaranteed to be headers but may not have equal contents
        //

        for (equivNo=sliceNo+1; equivNo<=maxNo; equivNo++)
        {
          sampleHeader = GPSmartVoxelProxy_GetHeader(This->fslices[equivNo]);
          if ( GPSmartVoxelHeader_operator_equal(sampleHeader,equivHeader) )
          {

            // Delete sampleHeader + proxy and replace with equivHeader/Proxy
            //
            GPSmartVoxelHeader_Destructor( sampleHeader );
            delete sampleHeader;
            delete This->fslices[equivNo];
            This->fslices[equivNo] = equivProxy;
          }
          else
          {
            // Not equal. Set this header to be
            // the current header for comparisons
            //
            equivProxy  = This->fslices[equivNo];
            equivHeader = GPSmartVoxelProxy_GetHeader(equivProxy);
          }

        }
        // Skip past examined slices
        //
        sliceNo = maxNo;
      }
    }
  }
}


// ***************************************************************************
// Calculate a "quality value" for the specified vector of voxels.
// The value returned should be >0 and such that the smaller the number
// the higher the quality of the slice.
//
// Preconditions: pSlice must consist of GPSmartVoxelNodeProxies only
// Process:
// o Examine each node in turn, summing:
//      no. of non-empty nodes
//      no. of volumes in each node
// o Calculate Quality=sigma(volumes in nod)/(no. of non-empty nodes)
//      if all nodes empty, return kInfinity
// o Call G4Exception on finding a GPSmartVoxelHeaderProxy
// ***************************************************************************
//
G4double GPSmartVoxelHeader_CalculateQuality(
		GPSmartVoxelHeader *This,
		GPProxyVector *pSlice)
{
  (void)This;
	
  G4double quality;
  G4int nNodes = pSlice->size();
  G4int noContained, maxContained=0, sumContained=0, sumNonEmptyNodes=0;
  GPSmartVoxelNode *node;

  for (G4int i=0; i<nNodes; i++)
  {
    if (GPSmartVoxelProxy_IsNode((*pSlice)[i]))
    {
      // Definitely a node. Add info to running totals
      //
      node = GPSmartVoxelProxy_GetNode((*pSlice)[i]);
      noContained = GPSmartVoxelNode_GetNoContained(node);
      if (noContained)
      {
        sumNonEmptyNodes++;
        sumContained += noContained;
        //
        // Calc maxContained for statistics
        //
        if (noContained>maxContained)
        {
          maxContained = noContained;
        }
      }
    }
  }

  // Calculate quality with protection against no non-empty nodes
  //
  if (sumNonEmptyNodes)
  {
    quality = sumContained/sumNonEmptyNodes;
  }
  else
  {
    quality = kInfinity;
  }

  return quality;
}

// ***************************************************************************
// Calculates and stores the minimum and maximum equivalent neighbour
// values for all slices at our level.
//
// Precondition: all slices are nodes.
// For each potential start of a group of equivalent nodes:
// o searches forwards in fslices to find group end
// o loops from start to end setting start and end slices.
// ***************************************************************************
//
void GPSmartVoxelHeader_BuildEquivalentSliceNos( GPSmartVoxelHeader *This )
{
  G4int sliceNo, minNo, maxNo, equivNo;
  G4int maxNode = GPSmartVoxelHeader_GetNoSlices(This);
  GPSmartVoxelNode *startNode, *sampleNode;
  for (sliceNo=0; sliceNo<maxNode; sliceNo++)
  {
    minNo = sliceNo;

    // Get first node (see preconditions - will throw exception if a header)
    //
    startNode = GPSmartVoxelProxy_GetNode(This->fslices[minNo]);
    assert(startNode != NULL);

    // Find max equivalent
    //
    for (equivNo=minNo+1; equivNo<maxNode; equivNo++)
    {
      sampleNode = GPSmartVoxelProxy_GetNode(This->fslices[equivNo]);
      assert( sampleNode != NULL );
      if (!(GPSmartVoxelNode_operator_equal(startNode,sampleNode))) { break; }
    }
    maxNo = equivNo-1;
    if (maxNo != minNo)
    {
      // Set min and max nos
      //
      for (equivNo=minNo; equivNo<=maxNo; equivNo++)
      {
        sampleNode = GPSmartVoxelProxy_GetNode(This->fslices[equivNo]);
        GPSmartVoxelNode_SetMinEquivalentSliceNo(sampleNode,minNo);
        GPSmartVoxelNode_SetMaxEquivalentSliceNo(sampleNode,maxNo);
      }
      // Advance outer loop to end of equivalent group
      //
      sliceNo = maxNo;
    }
  }
}

// ***************************************************************************
// Examined each contained node, refines (creates a replacement additional
// dimension of voxels) when there is more than one voxel in the slice.
// Does not refine further if already limited in two dimensions (=> this
// is the third level of limits)
//
// Preconditions: slices (nodes) have been built.
// ***************************************************************************
//
void GPSmartVoxelHeader_RefineNodes(
		GPSmartVoxelHeader *This,
		GPLogicalVolume* pVolume,
		GPVoxelLimits pLimits,
		G4double smartless)
{
  G4int refinedDepth=0, minVolumes;
  G4int maxNode = GPSmartVoxelHeader_GetNoSlices(This);

  if (GPVoxelLimits_IsXLimited(&pLimits))
  {
    refinedDepth++;
  }
  if (GPVoxelLimits_IsYLimited(&pLimits)) 
  {
    refinedDepth++;
  }
  if (GPVoxelLimits_IsZLimited(&pLimits))
  {
    refinedDepth++;
  }
  
  //const G4int kMinVoxelVolumesLevel2 = K_MIN_VOXEL_VOLUMES_LEVEL_2;
  //const G4int kMinVoxelVolumesLevel3 = K_MIN_VOXEL_VOLUMES_LEVEL_3;

  // Calculate minimum number of volumes necessary to refine
  //
  switch (refinedDepth)
  {
    case 0:
      minVolumes=minVoxelVolumesLevel2;
      break;
    case 1:
      minVolumes=minVoxelVolumesLevel3;
      break;
    default:
      minVolumes=10000;   // catch refinedDepth=3 and errors
      break;
  }

  if (refinedDepth<2)
  {
    G4int targetNo, noContainedDaughters, minNo, maxNo, replaceNo, i;
    G4double sliceWidth = (This->fmaxExtent-This->fminExtent)/maxNode;
    GPVoxelLimits newLimits; GPVoxelLimits_Constructor(&newLimits);
    GPSmartVoxelNode* targetNode;
    GPSmartVoxelProxy* targetNodeProxy;
    GPSmartVoxelHeader* replaceHeader;
    GPSmartVoxelProxy* replaceHeaderProxy;
    GPVolumeNosVector* targetList;
    GPSmartVoxelProxy* lastProxy;
      
    for (targetNo=0; targetNo<maxNode; targetNo++)
    {
      // Assume all slices are nodes (see preconditions)
      //
      targetNodeProxy = This->fslices[targetNo];
      targetNode = GPSmartVoxelProxy_GetNode(targetNodeProxy);

      if (GPSmartVoxelNode_GetNoContained(targetNode) >= minVolumes)
      {
        noContainedDaughters = GPSmartVoxelNode_GetNoContained(targetNode);
        targetList = new GPVolumeNosVector();
        targetList->reserve(noContainedDaughters);
        
        for (i=0; i<noContainedDaughters; i++)
        {
          targetList->push_back(GPSmartVoxelNode_GetVolume(targetNode,i));
        }
        minNo = GPSmartVoxelNode_GetMinEquivalentSliceNo(targetNode);
        maxNo = GPSmartVoxelNode_GetMaxEquivalentSliceNo(targetNode);

        // Delete node proxies at start of collected sets of nodes/headers
        //
        lastProxy=0;
        for (replaceNo=minNo; replaceNo<=maxNo; replaceNo++)
        {
          if (lastProxy != This->fslices[replaceNo])
          {
            lastProxy=This->fslices[replaceNo];
            delete lastProxy;
          }
        }
        // Delete node to be replaced
        //
        GPSmartVoxelNode_Destructor( targetNode );
        delete targetNode;

        // Create new headers + proxies and replace in fslices
        //
        newLimits = pLimits;
        GPVoxelLimits_AddLimit(&newLimits,This->faxis,This->fminExtent+sliceWidth*minNo,
                           This->fminExtent+sliceWidth*(maxNo+1));
                           
        replaceHeader = new GPSmartVoxelHeader;
        GPSmartVoxelHeader_Constructor2(
			replaceHeader,pVolume,newLimits,targetList,replaceNo,smartless);
			
        GPSmartVoxelHeader_SetMinEquivalentSliceNo(replaceHeader,minNo);
        GPSmartVoxelHeader_SetMaxEquivalentSliceNo(replaceHeader,maxNo);
        replaceHeaderProxy = new GPSmartVoxelProxy;
        GPSmartVoxelProxy_SetHeader(replaceHeaderProxy,replaceHeader);
        
        for (replaceNo=minNo; replaceNo<=maxNo; replaceNo++)
        {
          This->fslices[replaceNo] = replaceHeaderProxy;
        }
        // Finished replacing current `equivalent' group
        //
        delete targetList;
        targetNo=maxNo;
      }
    }
  }
}

// ***************************************************************************
// Builds and refines voxels between specified limits, considering only
// the physical volumes numbered `pCandidates'.
// o Chooses axis
// o Determines min and max extents (of mother solid) within limits.
// ***************************************************************************
//
void
GPSmartVoxelHeader_BuildVoxelsWithinLimits(GPSmartVoxelHeader *This, 
					   GPLogicalVolume* pVolume,
					   GPVoxelLimits pLimits,
					   const GPVolumeNosVector* pCandidates,
					   G4double smartless)
{
  // Choose best axis for slicing by:
  // 1. Trying all unlimited cartesian axes
  // 2. Select axis which gives greatest no slices
  
  GPProxyVector *pGoodSlices = NULL, *pTestSlices, *tmpSlices;
  G4double goodSliceScore=kInfinity, testSliceScore;
  EAxis goodSliceAxis = kXAxis;
  EAxis testAxis      = kXAxis;
  G4int node, maxNode, iaxis;
  GPVoxelLimits noLimits; GPVoxelLimits_Constructor(&noLimits);

  // Try all non-limited cartesian axes
  //
  for (iaxis=0; iaxis<3; iaxis++)
  {
    switch(iaxis)
    {
      case 0:
        testAxis = kXAxis;
        break;
      case 1:
        testAxis = kYAxis;
        break;
      case 2:
        testAxis = kZAxis;
        break;
    }
    if (!GPVoxelLimits_IsLimited2(&pLimits,testAxis))
    {
      pTestSlices = GPSmartVoxelHeader_BuildNodes(This,pVolume,pLimits,pCandidates,testAxis,smartless);
      
      if (pTestSlices == NULL) return;
      
      testSliceScore = GPSmartVoxelHeader_CalculateQuality(This,pTestSlices);
      if ( (!pGoodSlices) || (testSliceScore<goodSliceScore) )
      {
        goodSliceAxis  = testAxis;
        goodSliceScore = testSliceScore;
        tmpSlices      = pGoodSlices;
        pGoodSlices    = pTestSlices;
        pTestSlices    = tmpSlices;
      }
      if (pTestSlices)
      {
        // Destroy pTestSlices and all its contents
        //
        maxNode=pTestSlices->size();
        for (node=0; node<maxNode; node++)
        {
	  GPSmartVoxelNode *dyingNode = GPSmartVoxelProxy_GetNode((*pTestSlices)[node]);
	  GPSmartVoxelNode_Destructor( dyingNode );
	  delete dyingNode;
        }
        GPSmartVoxelProxy* tmpProx;
        while (pTestSlices->size()>0)
        {
          tmpProx = pTestSlices->back();
          pTestSlices->pop_back();
          for (GPProxyVector::iterator i=pTestSlices->begin();
                                       i!=pTestSlices->end(); i++)
          {
            if (*i==tmpProx)
            {
              pTestSlices->erase(i); i--;
            }
          }
          if ( tmpProx ) { delete tmpProx; }
        } 
        delete pTestSlices;
      }
    }
  }
  // Check for error case.. when limits already 3d,
  // so cannot select a new axis
  //
  assert( pGoodSlices != NULL );

  // 
  // We have selected pGoodSlices, with a score testSliceScore
  //

  // Store chosen axis, slice ptr
  //
  
  // TODO!!
  GPSmartVoxelHeader_SetSlices(This,pGoodSlices);
  
  delete pGoodSlices;   // Destroy slices vector, but not contained
                        // proxies or nodes
  This->faxis=goodSliceAxis;

  // Calculate and set min and max extents given our axis
  //
  
  GPVSolid* outerSolid = GPLogicalVolume_GetSolid(pVolume);
  const GPAffineTransform origin = GPAffineTransform_create_id();
  
  assert( outerSolid != NULL );
  
  if(!GPVSolid_CalculateExtent(outerSolid,This->faxis,pLimits,origin,&(This->fminExtent),&(This->fmaxExtent)))
  {
    GPVSolid_CalculateExtent(outerSolid,This->faxis,noLimits,origin,&(This->fminExtent),&(This->fmaxExtent));
  }

  // Calculate equivalent nos
  //
  GPSmartVoxelHeader_BuildEquivalentSliceNos(This);
  GPSmartVoxelHeader_CollectEquivalentNodes(This);     // Collect common nodes
  GPSmartVoxelHeader_RefineNodes(This,pVolume,pLimits,smartless); // Refine nodes creating headers

  // No common headers can exist because collapsed by construction
}

// ***************************************************************************
// Builds voxels for daughters specified volume, in NON-REPLICATED case
// o Create List of target volume nos (all daughters; 0->noDaughters-1)
// o BuildWithinLimits does Build & also determines mother dimensions.
// ***************************************************************************
//
void GPSmartVoxelHeader_BuildVoxels(GPSmartVoxelHeader *This, 
				    GPLogicalVolume* pVolume, 
				    G4double smartless)
{
  GPVoxelLimits limits;   // Create `unlimited' limits object
  GPVoxelLimits_Constructor(&limits);
  G4int nDaughters = GPLogicalVolume_GetNoDaughters(pVolume);
  
  GPVolumeNosVector targetList;
  targetList.reserve(nDaughters);
  
  for (G4int i=0; i<nDaughters; i++)
    {
      targetList.push_back(i);
    }
  GPSmartVoxelHeader_BuildVoxelsWithinLimits(This, pVolume, limits, &targetList, smartless);
}


// ***************************************************************************
// Returns true if all slices have equal contents.
// Preconditions: all equal slices have been collected.
// Procedure:
// o checks all slice proxy pointers are equal
// o returns true if only one slice or all slice proxies pointers equal.
// ***************************************************************************
//
G4bool GPSmartVoxelHeader_AllSlicesEqual( GPSmartVoxelHeader *This )
{
  G4int noSlices = GPSmartVoxelHeader_GetNoSlices(This);
  GPSmartVoxelProxy* refProxy;

  if (noSlices>1)
  {
    refProxy=This->fslices[0];
    for (G4int i=1; i<noSlices; i++)
    {
      if (refProxy!=This->fslices[i])
      {
        return false;
      }
    }
  }
  return true;
}

//}


void GPSmartVoxelNode_Constructor( GPSmartVoxelNode *This, G4int no )
{
  This->fmaxEquivalent = no;
  This->fminEquivalent = no;
  This->fcontents = NULL;
  This->fnumContents = 0;
}

void GPSmartVoxelNode_Destructor( GPSmartVoxelNode *This )
{
  std::free( This->fcontents );
}

void GPSmartVoxelNode_Insert( GPSmartVoxelNode *This, G4int thing )
{
  This->fnumContents++;
  This->fcontents = (G4int*)std::realloc( This->fcontents, sizeof(G4int)*This->fnumContents );
  This->fcontents[This->fnumContents-1] = thing;
}

void GPSmartVoxelProxy_SetNode( GPSmartVoxelProxy *This, GPSmartVoxelNode *n )
{
  This->fNode = n;
  This->fHeader = NULL;
}

void GPSmartVoxelProxy_SetHeader( GPSmartVoxelProxy *This, GPSmartVoxelHeader *h )
{
  This->fNode = NULL;
  This->fHeader = h;
}

#endif
