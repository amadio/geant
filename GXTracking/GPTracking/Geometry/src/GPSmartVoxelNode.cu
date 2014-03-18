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
// $Id: G4SmartVoxelNode.cc,v 1.6 2006-06-29 18:33:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class G4SmartVoxelNode
//
// Implementation
//
// --------------------------------------------------------------------

#include "GPSmartVoxelNode.h"

// Empty destructor
//
/*
FQUALIFIER
void GPSmartVoxelNode_Destructor()
{
  ;
}
*/

// Return true if contents equal
//
// Preconditions:
//
// Node contents were entered in the same order
FQUALIFIER
G4bool GPSmartVoxelNode_operator_equal (GPSmartVoxelNode *This,
					GPSmartVoxelNode *v)
{
  G4int maxNode=GPSmartVoxelNode_GetNoContained(This);
  if (maxNode==GPSmartVoxelNode_GetNoContained(v))
  {
    for (G4int node=0;node<maxNode;node++)
    {
      if (GPSmartVoxelNode_GetVolume(This,node) !=
	  GPSmartVoxelNode_GetVolume(v,node))
      {
        return false;
      }
    }
    return true;
  }
  return false;
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
// $Id: G4SmartVoxelNode.icc,v 1.7 2008-01-24 15:47:23 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4SmartVoxelNode Inline implementation
//
// --------------------------------------------------------------------

FQUALIFIER
G4int GPSmartVoxelNode_GetVolume(GEOMETRYLOC GPSmartVoxelNode *This,
				 G4int pVolumeNo)
{
  //check the validity of pVolumeNo
  return This->fcontents[pVolumeNo];
}

FQUALIFIER
G4int GPSmartVoxelNode_GetNoContained(GEOMETRYLOC GPSmartVoxelNode *This)
{
  return This->fnumContents;
}

/*
FQUALIFIER
void GPSmartVoxelNode_Insert(GPSmartVoxelNode *This,
			     G4int pVolumeNo)
{
  fcontents.push_back(pVolumeNo);
}

FQUALIFIER
G4int GPSmartVoxelNode_GetCapacity(GPSmartVoxelNode *This,) const
{
  return fcontents.capacity();
}

FQUALIFIER
void GPSmartVoxelNode_Reserve(GPSmartVoxelNode *This, G4int noSlices)
{
  fcontents.reserve(noSlices);
}

FQUALIFIER
void GPSmartVoxelNode_Shrink(GPSmartVoxelNode *This,)
{
  G4SliceVector(fcontents).swap(fcontents);
}
*/

FQUALIFIER
G4int GPSmartVoxelNode_GetMaxEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelNode *This)
{
  return This->fmaxEquivalent;
}

FQUALIFIER
void GPSmartVoxelNode_SetMaxEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelNode *This,
					      G4int pMax)
{
  This->fmaxEquivalent=pMax;
}

FQUALIFIER
G4int GPSmartVoxelNode_GetMinEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelNode *This)
{
  return This->fminEquivalent;
}

FQUALIFIER
void GPSmartVoxelNode_SetMinEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelNode *This,
					      G4int pMin)
{
  This->fminEquivalent=pMin;
}
