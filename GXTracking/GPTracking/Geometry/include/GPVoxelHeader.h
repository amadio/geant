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
// $Id: G4SmartVoxelHeader.hh,v 1.10 2006-06-29 18:32:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4SmartVoxelHeader
//
// Class description:
//
// Represents a set of voxels, created by a single axis of virtual division.
// Contains the individual voxels, which are potentially further divided
// along different axes.
//
// Member data:
//
// EAxis faxis
//   - The (cartesian) slicing/division axis
// G4double fmaxExtent
// G4double fminExtent
//   - Minimum and maximum coordiantes along the axis
// std::vector<G4SmartVoxelProxy*> fslices
//   - The slices along the axis
//
// G4int fminEquivalent
// G4int fmaxEquivalent
//   - Minimum and maximum equivalent slice nos.
//     [Applies to the level of the header, not its nodes]

// History:
// 18.04.01 G.Cosmo Migrated to STL vector
// 13.07.95 P.Kent  Initial version
// --------------------------------------------------------------------
#ifndef GPVOXELHEADER_HH
#define GPVOXELHEADER_HH

#include "GPVoxeldefs.h"
#include "GPVoxelLimits.h"
#include "GPLogicalVolume.h"
#include "GPVPhysicalVolume.h"
#include "GPAffineTransform.h"

#include "GPSmartVoxelHeader.h"
#include "GPSmartVoxelNode.h"
#include "GPSmartVoxelProxy.h"

namespace 
{
  typedef std::vector<GPSmartVoxelProxy*> GPProxyVector;
  typedef std::vector<GPSmartVoxelNode*> GPNodeVector;
  typedef std::vector<G4int> GPVolumeNosVector;
  typedef std::vector<G4double> GPVolumeExtentVector;
}

extern "C" {

static int minVoxelVolumesLevel1 = kMinVoxelVolumesLevel1;
static int minVoxelVolumesLevel2 = kMinVoxelVolumesLevel2;
static int minVoxelVolumesLevel3 = kMinVoxelVolumesLevel3;

void GPSmartVoxelHeader_SetMinVoxelLimits( int lev1, int lev2, int lev3 );
int GPSmartVoxelHeader_GetMinVoxelVolumesLevel1();
int GPSmartVoxelHeader_GetMinVoxelVolumesLevel2();
int GPSmartVoxelHeader_GetMinVoxelVolumesLevel3();

void GPSmartVoxelNode_Destructor( GPSmartVoxelNode *This );
void GPSmartVoxelNode_Constructor( GPSmartVoxelNode *This, G4int no );
void GPSmartVoxelNode_Insert( GPSmartVoxelNode *This, G4int thing );
void GPSmartVoxelProxy_SetNode( GPSmartVoxelProxy *This, GPSmartVoxelNode *n );
void GPSmartVoxelProxy_SetHeader( GPSmartVoxelProxy *This, GPSmartVoxelHeader *h );

void GPSmartVoxelHeader_Constructor2(GPSmartVoxelHeader *This,
                                     GPLogicalVolume* pVolume,
                                     GPVoxelLimits pLimits,
                                     const GPVolumeNosVector* pCandidates,
                                     G4int pSlice,
                                     G4double smartless);

void GPSmartVoxelHeader_BuildVoxels(GPSmartVoxelHeader *This, 
                                    GPLogicalVolume* pVolume, 
                                    G4double smartless);
void
GPSmartVoxelHeader_BuildVoxelsWithinLimits(GPSmartVoxelHeader *This, 
                                           GPLogicalVolume* pVolume,
                                           GPVoxelLimits pLimits,
                                           const GPVolumeNosVector* pCandidates,
                                           G4double smartless);


void GPSmartVoxelHeader_Constructor(GPSmartVoxelHeader *This,
                                    GPLogicalVolume* pVolume,
                                    G4int pSlice,
                                    G4double smartless );

void GPSmartVoxelHeader_Constructor2(GPSmartVoxelHeader *This,
                                     GPLogicalVolume* pVolume,
                                     GPVoxelLimits pLimits,
                                     const GPVolumeNosVector* pCandidates,
                                     G4int pSlice,
                                     G4double smartless);

void
GPSmartVoxelHeader_Destructor(GPSmartVoxelHeader *This);

G4bool GPSmartVoxelHeader_operator_equal(GPSmartVoxelHeader *This,
                                         GPSmartVoxelHeader* pHead);

void GPSmartVoxelHeader_SetSlices(GPSmartVoxelHeader *This,
                                  const GPProxyVector *slices );

GPProxyVector* GPSmartVoxelHeader_BuildNodes(GPSmartVoxelHeader *This,
                                             GPLogicalVolume* pVolume,
                                             GPVoxelLimits pLimits,
                                             const GPVolumeNosVector* pCandidates,
                                             EAxis pAxis,
                                             G4double smartlessUser );

void GPSmartVoxelHeader_CollectEquivalentNodes( GPSmartVoxelHeader *This );

void GPSmartVoxelHeader_CollectEquivalentHeaders( GPSmartVoxelHeader *This );

G4double GPSmartVoxelHeader_CalculateQuality(GPSmartVoxelHeader *This,
					     GPProxyVector *pSlice);

void GPSmartVoxelHeader_BuildEquivalentSliceNos( GPSmartVoxelHeader *This );

void GPSmartVoxelHeader_RefineNodes(GPSmartVoxelHeader *This,
				    GPLogicalVolume* pVolume,
				    GPVoxelLimits pLimits,
				    G4double smartless);

void
GPSmartVoxelHeader_BuildVoxelsWithinLimits(GPSmartVoxelHeader *This, 
                                           GPLogicalVolume* pVolume,
                                           GPVoxelLimits pLimits,
                                           const GPVolumeNosVector* pCandidates,
                                           G4double smartless);

void GPSmartVoxelHeader_BuildVoxels(GPSmartVoxelHeader *This, 
                                    GPLogicalVolume* pVolume, 
                                    G4double smartless);

G4bool GPSmartVoxelHeader_AllSlicesEqual( GPSmartVoxelHeader *This );

void GPSmartVoxelNode_Constructor( GPSmartVoxelNode *This, G4int no );

void GPSmartVoxelNode_Destructor( GPSmartVoxelNode *This );

void GPSmartVoxelNode_Insert( GPSmartVoxelNode *This, G4int thing );

void GPSmartVoxelProxy_SetNode( GPSmartVoxelProxy *This, GPSmartVoxelNode *n );

void GPSmartVoxelProxy_SetHeader( GPSmartVoxelProxy *This, GPSmartVoxelHeader *h );

}

#endif
