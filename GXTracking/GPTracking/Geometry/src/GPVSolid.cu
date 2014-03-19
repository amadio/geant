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
// $Id: GPVSolid.cc,v 1.40 2010-10-19 15:19:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class GPVSolid
//
// Implementation for solid base class
//
// History:
//
//  06.12.02 V.Grichine, restored original conditions in ClipPolygon()
//  10.05.02 V.Grichine, ClipPolygon(): clip only other axis and limited voxels
//  15.04.02 V.Grichine, bug fixed in ClipPolygon(): clip only one axis
//  13.03.02 V.Grichine, cosmetics of voxel limit functions  
//  15.11.00 D.Williams, V.Grichine, fix in CalculateClippedPolygonExtent()
//  10.07.95 P.Kent, Added == operator, solid Store entry
//  30.06.95 P.Kent, Created.
// --------------------------------------------------------------------

#include "GPVSolid.h"
//#include "G4SolidStore.hh"
//#include "globals.hh"
//#include "Randomize.hh"
//#include "G4GeometryTolerance.hh"

#include "GPVoxelLimits.h"
#include "GPAffineTransform.h"
//#include "G4VisExtent.hh"

#include "GPBox.h"
#include "GPTrd.h"
#include "GPTubs.h"
#include "GPCons.h"
#include "GPOrb.h"

//////////////////////////////////////////////////////////////////////////
//
// Constructor
//  - Copies name
//  - Add ourselves to solid Store

FQUALIFIER
void GPVSolid_Constructor(GPVSolid *This,
                          ESolid kType)
{
  This->fType = kType;
  // Register to store
  //
  //    G4SolidStore::GetInstance()->Register(this);
}

//////////////////////////////////////////////////////////////////////////
//
// Dummy implementations ...
/*
const GPVSolid* GPVSolid::GetConstituentSolid(G4int) const
{ return 0; } 

GPVSolid* GPVSolid::GetConstituentSolid(G4int)
{ return 0; } 

const G4DisplacedSolid* GPVSolid::GetDisplacedSolidPtr() const
{ return 0; } 

G4DisplacedSolid* GPVSolid::GetDisplacedSolidPtr() 
{ return 0; } 
*/
////////////////////////////////////////////////////////////////
//
// Returns an estimation of the solid volume in internal units.
// The number of statistics and error accuracy is fixed.
// This method may be overloaded by derived classes to compute the
// exact geometrical quantity for solids where this is possible.
// or anyway to cache the computed value.
// This implementation does NOT cache the computed value.

///////////////////////////////////////////////////////////////////////////
// 
// Calculate the maximum and minimum extents of the polygon described
// by the vertices: pSectionIndex->pSectionIndex+1->
//                   pSectionIndex+2->pSectionIndex+3->pSectionIndex
// in the List pVertices
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices
//

FQUALIFIER
void GPVSolid_ClipCrossSection( GPThreeVectorList* pVertices,
                                const G4int pSectionIndex,
                                GPVoxelLimits* pVoxelLimit,
                                const EAxis pAxis, 
				G4double *pMin, G4double *pMax)
{
  GPThreeVectorList polygon;
  GPThreeVectorList_Constructor(&polygon);
  //  polygon.reserve(4);
  //  polygon.push_back((*pVertices)[pSectionIndex]);
  //  polygon.push_back((*pVertices)[pSectionIndex+1]);
  //  polygon.push_back((*pVertices)[pSectionIndex+2]);
  //  polygon.push_back((*pVertices)[pSectionIndex+3]);
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+1));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+2));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+3));

  //  G4cout<<"ClipCrossSection: 0-1-2-3"<<G4endl;
  GPVSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  return;
}

//////////////////////////////////////////////////////////////////////////////////
//
// Calculate the maximum and minimum extents of the polygons
// joining the CrossSections at pSectionIndex->pSectionIndex+3 and
//                              pSectionIndex+4->pSectionIndex7
//
// in the List pVertices, within the boundaries of the voxel limits pVoxelLimit
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices

FQUALIFIER
void GPVSolid_ClipBetweenSections( GPThreeVectorList* pVertices,
                                   const G4int pSectionIndex,
                                   GPVoxelLimits* pVoxelLimit,
                                   const EAxis pAxis, 
				   G4double *pMin, G4double *pMax)
{
  GPThreeVectorList polygon;
  GPThreeVectorList_Constructor(&polygon);
  /*
  polygon.reserve(4);
  polygon.push_back((*pVertices)[pSectionIndex]);
  polygon.push_back((*pVertices)[pSectionIndex+4]);
  polygon.push_back((*pVertices)[pSectionIndex+5]);
  polygon.push_back((*pVertices)[pSectionIndex+1]);
  */
  // G4cout<<"ClipBetweenSections: 0-4-5-1"<<G4endl;

  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+4));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+5));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+1));

  GPVSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);

  /*
  polygon.clear();
  polygon.push_back((*pVertices)[pSectionIndex+1]);
  polygon.push_back((*pVertices)[pSectionIndex+5]);
  polygon.push_back((*pVertices)[pSectionIndex+6]);
  polygon.push_back((*pVertices)[pSectionIndex+2]);
  */
  // G4cout<<"ClipBetweenSections: 1-5-6-2"<<G4endl;

  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+1));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+5));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+6));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+2));

  GPVSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  /*
  polygon.clear();
  polygon.push_back((*pVertices)[pSectionIndex+2]);
  polygon.push_back((*pVertices)[pSectionIndex+6]);
  polygon.push_back((*pVertices)[pSectionIndex+7]);
  polygon.push_back((*pVertices)[pSectionIndex+3]);
  */
  //  G4cout<<"ClipBetweenSections: 2-6-7-3"<<G4endl;

  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+2));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+6));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+7));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+3));

  GPVSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  /*
  polygon.clear();
  polygon.push_back((*pVertices)[pSectionIndex+3]);
  polygon.push_back((*pVertices)[pSectionIndex+7]);
  polygon.push_back((*pVertices)[pSectionIndex+4]);
  polygon.push_back((*pVertices)[pSectionIndex]);
  */
  //  G4cout<<"ClipBetweenSections: 3-7-4-0"<<G4endl;

  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+3));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+7));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex+4));
  GPThreeVectorList_pushback(&polygon,GPThreeVectorList_popback(pVertices,pSectionIndex));

  GPVSolid_CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  return;
}


///////////////////////////////////////////////////////////////////////////////
//
// Calculate the maximum and minimum extents of the convex polygon pPolygon
// along the axis pAxis, within the limits pVoxelLimit
//

FQUALIFIER
void
GPVSolid_CalculateClippedPolygonExtent(GPThreeVectorList pPolygon,
				       GPVoxelLimits* pVoxelLimit,
				       const EAxis pAxis, 
				       G4double *pMin,
				       G4double *pMax) 
{
  G4int noLeft,i;
  G4double component;
  /*  
  G4cout<<G4endl;
  for(i = 0 ; i < pPolygon.size() ; i++ )
  {
      G4cout << i << "\t"
             << "p.x = " << pPolygon[i].operator()(pAxis) << "\t"
        //   << "p.y = " << pPolygon[i].y() << "\t"
        //   << "p.z = " << pPolygon[i].z() << "\t"
             << G4endl;
  }    
  G4cout<<G4endl;
  */  
  GPVSolid_ClipPolygon(pPolygon,pVoxelLimit,pAxis);
  noLeft = GPThreeVectorList_size(&pPolygon);

  if ( noLeft )
  {
    //  G4cout<<G4endl;
    for (i=0;i<noLeft;i++)
    {
      //      component = pPolygon[i].operator()(pAxis);
      switch(pAxis) {
      case kXAxis: 
	//	component = pPolygon[i].x;
	component = (GPThreeVectorList_popback(&pPolygon,i)).x;
      case kYAxis: 
	//	component = pPolygon[i].y;
	component = (GPThreeVectorList_popback(&pPolygon,i)).y;
      case kZAxis: 
	//	component = pPolygon[i].z;
	component = (GPThreeVectorList_popback(&pPolygon,i)).z;
      default: 
	component = 0.0;
      }
      //  G4cout <<i<<"\t"<<component<<G4endl;
 
      if (component < *pMin) 
      { 
        //  G4cout <<i<<"\t"<<"Pmin = "<<component<<G4endl;
        *pMin = component;      
      }
      if (component > *pMax)
      {  
        //  G4cout <<i<<"\t"<<"PMax = "<<component<<G4endl;
        *pMax = component;  
      }    
    }
    //  G4cout<<G4endl;
  }
  // G4cout<<"pMin = "<<pMin<<"\t"<<"pMax = "<<pMax<<G4endl;
}

/////////////////////////////////////////////////////////////////////////////
//
// Clip the convex polygon described by the vertices at
// pSectionIndex ->pSectionIndex+3 within pVertices to the limits pVoxelLimit
//
// Set pMin to the smallest
//
// Calculate the extent of the polygon along pAxis, when clipped to the
// limits pVoxelLimit. If the polygon exists after clippin, set pMin to
// the polygon's minimum extent along the axis if <pMin, and set pMax to
// the polygon's maximum extent along the axis if >pMax.
//
// The polygon is described by a set of vectors, where each vector represents
// a vertex, so that the polygon is described by the vertex sequence:
//   0th->1st 1st->2nd 2nd->... nth->0th
//
// Modifications to the polygon are made
//
// NOTE: Execessive copying during clipping

FQUALIFIER
void GPVSolid_ClipPolygon( GPThreeVectorList pPolygon,
                           GPVoxelLimits* pVoxelLimit,
                           const EAxis                        )
{
  GPThreeVectorList outputPolygon;
  GPThreeVectorList_Constructor(&outputPolygon);

  if ( GPVoxelLimits_IsLimited(pVoxelLimit) )
  {
    if ( GPVoxelLimits_IsXLimited(pVoxelLimit) ) // && pAxis != kXAxis)
    {
      GPVoxelLimits simpleLimit1;
      GPVoxelLimits_Constructor(&simpleLimit1);
      GPVoxelLimits_AddLimit(&simpleLimit1,kXAxis,
			     GPVoxelLimits_GetMinXExtent(pVoxelLimit),
			     kInfinity);
      //  G4cout<<"MinXExtent()"<<G4endl;
      GPVSolid_ClipPolygonToSimpleLimits(pPolygon,&outputPolygon,&simpleLimit1);
   
      //      pPolygon.clear();

      if ( ! GPThreeVectorList_size(&outputPolygon) )  return;

      GPVoxelLimits simpleLimit2;
      GPVoxelLimits_Constructor(&simpleLimit2);

      //  G4cout<<"MaxXExtent()"<<G4endl;
      GPVoxelLimits_AddLimit(&simpleLimit2,kXAxis,-kInfinity,
			     GPVoxelLimits_GetMaxXExtent(pVoxelLimit));
      GPVSolid_ClipPolygonToSimpleLimits(outputPolygon,&pPolygon,&simpleLimit2);

      if ( ! GPThreeVectorList_size(&pPolygon) )       return;
      else ; // outputPolygon.clear();
    }
    if ( GPVoxelLimits_IsYLimited(pVoxelLimit) ) // && pAxis != kYAxis)
    {
      GPVoxelLimits simpleLimit1;
      GPVoxelLimits_Constructor(&simpleLimit1);
      GPVoxelLimits_AddLimit(&simpleLimit1,kYAxis,
			     GPVoxelLimits_GetMinYExtent(pVoxelLimit),
			     kInfinity);
      GPVSolid_ClipPolygonToSimpleLimits(pPolygon,&outputPolygon,&simpleLimit1);

      // Must always clear pPolygon - for clip to simpleLimit2 and in case of

      //      pPolygon.clear();

      if ( ! GPThreeVectorList_size(&outputPolygon) )  return;

      GPVoxelLimits simpleLimit2;
      GPVoxelLimits_Constructor(&simpleLimit2);

      GPVoxelLimits_AddLimit(&simpleLimit2,kYAxis,-kInfinity,
			     GPVoxelLimits_GetMaxYExtent(pVoxelLimit));
      GPVSolid_ClipPolygonToSimpleLimits(outputPolygon,&pPolygon,&simpleLimit2);

      if ( ! GPThreeVectorList_size(&pPolygon) )       return;
      else ; //                        outputPolygon.clear();
    }
    if ( GPVoxelLimits_IsZLimited(pVoxelLimit) ) // && pAxis != kZAxis)
    {
      GPVoxelLimits simpleLimit1;
      GPVoxelLimits_Constructor(&simpleLimit1);

      GPVoxelLimits_AddLimit(&simpleLimit1,kZAxis,
			     GPVoxelLimits_GetMinZExtent(pVoxelLimit),
			     kInfinity);
      GPVSolid_ClipPolygonToSimpleLimits(pPolygon,&outputPolygon,&simpleLimit1);

      // Must always clear pPolygon - for clip to simpleLimit2 and in case of
      // early exit

      //      pPolygon.clear();

      if ( ! GPThreeVectorList_size(&outputPolygon) )  return;

      GPVoxelLimits simpleLimit2;
      GPVoxelLimits_Constructor(&simpleLimit2);

      GPVoxelLimits_AddLimit(&simpleLimit2,kZAxis,-kInfinity,
			     GPVoxelLimits_GetMaxZExtent(pVoxelLimit));
      GPVSolid_ClipPolygonToSimpleLimits(outputPolygon,&pPolygon,&simpleLimit2);

      // Return after final clip - no cleanup
    }
  }
}

////////////////////////////////////////////////////////////////////////////
//
// pVoxelLimits must be only limited along one axis, and either the maximum
// along the axis must be +kInfinity, or the minimum -kInfinity

FQUALIFIER
void
GPVSolid_ClipPolygonToSimpleLimits( GPThreeVectorList pPolygon,
				    GPThreeVectorList* outputPolygon,
				    GPVoxelLimits* pVoxelLimit       )
{
  G4int i;
  G4int noVertices=GPThreeVectorList_size(&pPolygon);
  GPThreeVector vEnd,vStart;

  for (i = 0 ; i < noVertices ; i++ )
  {
    vStart = GPThreeVectorList_popback(&pPolygon,i);
    // G4cout << "i = " << i << G4endl;
    if ( i == noVertices-1 )    vEnd = GPThreeVectorList_popback(&pPolygon,0);//  pPolygon[0];
    else                        vEnd = GPThreeVectorList_popback(&pPolygon,i+1);//pPolygon[i+1];

    if ( GPVoxelLimits_Inside(pVoxelLimit,vStart) )
    {
      if (GPVoxelLimits_Inside(pVoxelLimit,vEnd))
      {
        // vStart and vEnd inside -> output end point
        //
        GPThreeVectorList_pushback(outputPolygon,vEnd);
      }
      else
      {
        // vStart inside, vEnd outside -> output crossing point
        //
        // G4cout << "vStart inside, vEnd outside" << G4endl;
        GPVoxelLimits_ClipToLimits(pVoxelLimit,&vStart,&vEnd);
	//        outputPolygon.push_back(vEnd);
        GPThreeVectorList_pushback(outputPolygon,vEnd);
      }    
    }
    else
    {
      if (GPVoxelLimits_Inside(pVoxelLimit,vEnd))
      {
        // vStart outside, vEnd inside -> output inside section
        //
        // G4cout << "vStart outside, vEnd inside" << G4endl;
        GPVoxelLimits_ClipToLimits(pVoxelLimit,&vStart,&vEnd);
	//        outputPolygon.push_back(vStart);
	//        outputPolygon.push_back(vEnd);  
        GPThreeVectorList_pushback(outputPolygon,vStart);
        GPThreeVectorList_pushback(outputPolygon,vEnd);

      }
      else  // Both point outside -> no output
      {
        // outputPolygon.push_back(vStart);
        // outputPolygon.push_back(vEnd);  
      }
    }
  }
}

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
// $Id: GPVSolid.cc,v 1.40 2010-10-19 15:19:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class GPVSolid
//
// Implementation for solid base class
//
// History:
//
//  06.12.02 V.Grichine, restored original conditions in ClipPolygon()
//  10.05.02 V.Grichine, ClipPolygon(): clip only other axis and limited voxels
//  15.04.02 V.Grichine, bug fixed in ClipPolygon(): clip only one axis
//  13.03.02 V.Grichine, cosmetics of voxel limit functions  
//  15.11.00 D.Williams, V.Grichine, fix in CalculateClippedPolygonExtent()
//  10.07.95 P.Kent, Added == operator, solid Store entry
//  30.06.95 P.Kent, Created.
// --------------------------------------------------------------------

// Virtual Solid methods                                

FQUALIFIER
EInside GPVSolid_Inside( GEOMETRYLOC const GPVSolid *This, 
			 GPThreeVector p)
{
  switch(This->fType) {
  case kBox:
    return GPBox_Inside((GEOMETRYLOC GPBox*)This,p);
  case kTrd:
    return GPTrd_Inside((GEOMETRYLOC GPTrd*)This,p);
  case kTubs:
    return GPTubs_Inside((GEOMETRYLOC GPTubs*)This,p);
  case kCons:
    return GPCons_Inside((GEOMETRYLOC GPCons*)This,p);
  case kOrb:
    return GPOrb_Inside((GEOMETRYLOC GPOrb*)This,p);
  default:
    return kOutside;
  }
}

FQUALIFIER
GPThreeVector GPVSolid_SurfaceNormal(GEOMETRYLOC const GPVSolid *This, 
				     GPThreeVector p)
{
  switch(This->fType) {
  case kBox:
    return GPBox_SurfaceNormal((GEOMETRYLOC GPBox*)This,p);
  case kTrd:
    return GPTrd_SurfaceNormal((GEOMETRYLOC GPTrd*)This,p);
  case kTubs:
    return GPTubs_SurfaceNormal((GEOMETRYLOC GPTubs*)This,p);
  case kCons:
    return GPCons_SurfaceNormal((GEOMETRYLOC GPCons*)This,p);
  case kOrb:
    return GPOrb_SurfaceNormal((GEOMETRYLOC GPOrb*)This,p);
  default:
    return GPThreeVector_create(0,0,0);
  }
}

FQUALIFIER
G4double GPVSolid_DistanceToIn(GEOMETRYLOC const GPVSolid *This, 
			       GPThreeVector p)
{
  switch(This->fType) {
  case kBox:
    return GPBox_DistanceToIn((GEOMETRYLOC GPBox*)This,p);
  case kTrd:
    return GPTrd_DistanceToIn((GEOMETRYLOC GPTrd*)This,p);
  case kTubs:
    return GPTubs_DistanceToIn((GEOMETRYLOC GPTubs*)This,p);
  case kCons:
    return GPCons_DistanceToIn((GEOMETRYLOC GPCons*)This,p);
  case kOrb:
    return GPOrb_DistanceToIn((GEOMETRYLOC GPOrb*)This,p);
  default:
    return 0;

  }
}

FQUALIFIER
G4double GPVSolid_DistanceToOut(GEOMETRYLOC const GPVSolid *This, 
				GPThreeVector p)
{
  switch(This->fType) {
  case kBox:
    return GPBox_DistanceToOut((GEOMETRYLOC GPBox*)This,p);
  case kTrd:
    return GPTrd_DistanceToOut((GEOMETRYLOC GPTrd*)This,p);
  case kTubs:
    return GPTubs_DistanceToOut((GEOMETRYLOC GPTubs*)This,p);
  case kCons:
    return GPCons_DistanceToOut((GEOMETRYLOC GPCons*)This,p);
  case kOrb:
    return GPOrb_DistanceToOut((GEOMETRYLOC GPOrb*)This,p);
  default:
    return 0;
  }
}

FQUALIFIER
G4double GPVSolid_DistanceToIn2(GEOMETRYLOC const GPVSolid *This, 
				GPThreeVector p,
				GPThreeVector v)
{
  switch(This->fType) {
  case kBox:
    return GPBox_DistanceToIn2((GEOMETRYLOC GPBox*)This,p,v);
  case kTrd:
    return GPTrd_DistanceToIn2((GEOMETRYLOC GPTrd*)This,p,v);
  case kTubs:
    return GPTubs_DistanceToIn2((GEOMETRYLOC GPTubs*)This,p,v);
  case kCons:
    return GPCons_DistanceToIn2((GEOMETRYLOC GPCons*)This,p,v);
  case kOrb:
    return GPOrb_DistanceToIn2((GEOMETRYLOC GPOrb*)This,p,v);
  default:
    return 0;
  }
}

FQUALIFIER
G4double GPVSolid_DistanceToOut2(GEOMETRYLOC const GPVSolid *This,
				 GPThreeVector p,
				 GPThreeVector v,
				 const G4bool calcNorm,
				 G4bool *validNorm,
				 GPThreeVector *n)
{
  switch(This->fType) {
  case kBox:
    return GPBox_DistanceToOut2((GEOMETRYLOC GPBox*)This,p,v,calcNorm,validNorm,n);
  case kTrd:
    return GPTrd_DistanceToOut2((GEOMETRYLOC GPTrd*)This,p,v,calcNorm,validNorm,n);
  case kTubs:
    return GPTubs_DistanceToOut2((GEOMETRYLOC GPTubs*)This,p,v,calcNorm,validNorm,n);
  case kCons:
    return GPCons_DistanceToOut2((GEOMETRYLOC GPCons*)This,p,v,calcNorm,validNorm,n);
  case kOrb:
    return GPOrb_DistanceToOut2((GEOMETRYLOC GPOrb*)This,p,v,calcNorm,validNorm,n);
  default:
    return 0;
  }
}

FQUALIFIER
G4bool GPVSolid_CalculateExtent(const GPVSolid *This,
				const EAxis pAxis,
				GPVoxelLimits pVoxelLimit,
				GPAffineTransform pTransform,
				G4double* pMin, G4double* pMax)
{
  switch(This->fType) {
  case kBox:
    return GPBox_CalculateExtent((GEOMETRYLOC GPBox*)This, pAxis,
				 pVoxelLimit, pTransform, pMin, pMax);
  case kTrd:
    return GPTrd_CalculateExtent((GEOMETRYLOC GPTrd*)This, pAxis,
				  pVoxelLimit, pTransform, pMin, pMax);
  case kTubs:
    return GPTubs_CalculateExtent((GEOMETRYLOC GPTubs*)This, pAxis,
				  pVoxelLimit, pTransform, pMin, pMax);
  case kCons:
    return GPCons_CalculateExtent((GEOMETRYLOC GPCons*)This, pAxis,
				  pVoxelLimit, pTransform, pMin, pMax);
  case kOrb:
    return GPOrb_CalculateExtent((GEOMETRYLOC GPOrb*)This, pAxis,
				  pVoxelLimit, pTransform, pMin, pMax);
  default:
    return 0;
  }
}


