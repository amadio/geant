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
#ifndef GPSMARTVOXELHEADER_HH
#define GPSMARTVOXELHEADER_HH

#include "GPTypeDef.h"
#include "GPGeomdefs.h"

#include "GPSmartVoxelProxy.h"
#include "GPLogicalVolume.h"
#include "GPVoxelLimits.h"

//#include <vector>
// Typedefs
//typedef std::vector<G4SmartVoxelProxy*> G4ProxyVector;
//typedef std::vector<G4SmartVoxelNode*> G4NodeVector;
//typedef std::vector<G4int> G4VolumeNosVector;
//typedef std::vector<G4double> G4VolumeExtentVector;

struct GPSmartVoxelHeader
{
  G4int fminEquivalent;
  G4int fmaxEquivalent;
  // Min and max equivalent slice nos for previous level.

  EAxis faxis, fparamAxis;
  // Axis for slices.

  G4double fmaxExtent;
  G4double fminExtent;
  // Max and min coordinate along faxis.

  //  G4ProxyVector fslices;
  GEOMETRYLOC struct GPSmartVoxelProxy* *fslices;
  // Slices along axis.
  G4int fnumSlices;
};

extern "C" {

FQUALIFIER
G4int GPSmartVoxelHeader_GetMaxEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelHeader *This);

FQUALIFIER
void GPSmartVoxelHeader_SetMaxEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelHeader *This,
						G4int pMax);

FQUALIFIER
G4int GPSmartVoxelHeader_GetMinEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelHeader *This);

FQUALIFIER
void GPSmartVoxelHeader_SetMinEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelHeader *This,
						G4int pMin);

FQUALIFIER
EAxis GPSmartVoxelHeader_GetAxis(GEOMETRYLOC GPSmartVoxelHeader *This);

FQUALIFIER
EAxis GPSmartVoxelHeader_GetParamAxis(GEOMETRYLOC GPSmartVoxelHeader *This);

FQUALIFIER
G4double GPSmartVoxelHeader_GetMaxExtent(GEOMETRYLOC GPSmartVoxelHeader *This);

FQUALIFIER
G4double GPSmartVoxelHeader_GetMinExtent(GEOMETRYLOC GPSmartVoxelHeader *This);

FQUALIFIER
G4int GPSmartVoxelHeader_GetNoSlices(GEOMETRYLOC GPSmartVoxelHeader *This);

FQUALIFIER
GPSmartVoxelProxy* GPSmartVoxelHeader_GetSlice(GEOMETRYLOC GPSmartVoxelHeader *This,
					       G4int n);

}

#endif
