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
// $Id: G4SmartVoxelNode.hh,v 1.12 2008-01-24 15:47:23 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4SmartVoxelNode
//
// Class description:
//
// A node in the smart voxel hierarchy - a `slice' of space along a given
// axis between given minima and maxima. Note that the node is not aware
// of its position - this information being available/derivable by the
// node's owner(s) (voxelheaders).
//
// Member Data:
//
// G4int fminEquivalent
// G4int fmaxEquivalent
//   - Min and maximum nodes with same contents. Set by constructor
//     and set methods.
// std::vector<G4int> fcontents
//   - Vector of no.s of volumes inside the node.

// History:
// 18.04.01 G.Cosmo Migrated to STL vector
// 12.07.95 P.Kent  Initial version
// --------------------------------------------------------------------
#ifndef GPSMARTVOXELNODE_HH
#define GPSMARTVOXELNODE_HH

#include "GPTypeDef.h"

//#include <vector>
//typedef std::vector<G4int> G4SliceVector;

struct GPSmartVoxelNode
{
    G4int fminEquivalent;
    G4int fmaxEquivalent;
    //    G4SliceVector fcontents;
    GEOMETRYLOC G4int *fcontents;
    G4int fnumContents;
};

extern "C" {

FQUALIFIER
G4bool GPSmartVoxelNode_operator_equal (GPSmartVoxelNode *This,
					GPSmartVoxelNode *v);

FQUALIFIER
G4int GPSmartVoxelNode_GetVolume(GEOMETRYLOC GPSmartVoxelNode *This,
				 G4int pVolumeNo);

FQUALIFIER
G4int GPSmartVoxelNode_GetNoContained(GEOMETRYLOC GPSmartVoxelNode *This);

FQUALIFIER
G4int GPSmartVoxelNode_GetMaxEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelNode *This);

FQUALIFIER
void GPSmartVoxelNode_SetMaxEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelNode *This,
					      G4int pMax);

FQUALIFIER
G4int GPSmartVoxelNode_GetMinEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelNode *This);

FQUALIFIER
void GPSmartVoxelNode_SetMinEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelNode *This,
					      G4int pMin);
}
#endif
