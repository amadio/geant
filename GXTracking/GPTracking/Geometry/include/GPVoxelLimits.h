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
// $Id: G4VoxelLimits.hh,v 1.9 2006-06-29 18:33:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VoxelLimits
//
// Class description:
//
// Represents limitation/restrictions of space, where restrictions
// are only made perpendicular to the cartesian axes.
//
//
// Member data:
//
// G4double fxAxisMin,fxAxisMax
// G4double fyAxisMin,fyAxisMax
// G4double fzAxisMin,fzAxisMax
//   - The min and max values along each axis. +-kInfinity if not restricted.
//
//
// Notes:
//
// Beware no break statements after returns in switch(pAxis)s.

// History:
// 13.07.95 P.Kent Initial version.
// --------------------------------------------------------------------
#ifndef GPVOXELLIMITS_HH
#define GPVOXELLIMITS_HH

#include "GPTypeDef.h"
#include "GPGeomdefs.h"
#include "GPThreeVector.h"

struct GPVoxelLimits
{
  G4double fxAxisMin;
  G4double fxAxisMax;
  G4double fyAxisMin;
  G4double fyAxisMax;
  G4double fzAxisMin;
  G4double fzAxisMax;
};

extern "C" {

FQUALIFIER
void GPVoxelLimits_Constructor(GPVoxelLimits *This);

FQUALIFIER
void GPVoxelLimits_AddLimit( GPVoxelLimits *This,
                             const EAxis pAxis, 
                             const G4double pMin,
                             const G4double pMax );

FQUALIFIER
G4bool GPVoxelLimits_ClipToLimits( GPVoxelLimits *This, 
                                   GPThreeVector* pStart,
                                   GPThreeVector* pEnd      );

FQUALIFIER
G4int GPVoxelLimits_OutCode( GPVoxelLimits *This,
                             GPThreeVector pVec );

FQUALIFIER
G4double GPVoxelLimits_GetMaxXExtent(GPVoxelLimits *This );

FQUALIFIER
G4double GPVoxelLimits_GetMaxYExtent(GPVoxelLimits *This );

FQUALIFIER
G4double GPVoxelLimits_GetMaxZExtent(GPVoxelLimits *This );

FQUALIFIER
G4double GPVoxelLimits_GetMinXExtent(GPVoxelLimits *This );

FQUALIFIER
G4double GPVoxelLimits_GetMinYExtent(GPVoxelLimits *This );

FQUALIFIER
G4double GPVoxelLimits_GetMinZExtent(GPVoxelLimits *This );

FQUALIFIER
G4double GPVoxelLimits_GetMaxExtent(GPVoxelLimits *This ,
				    const EAxis pAxis);

FQUALIFIER
G4double GPVoxelLimits_GetMinExtent( GPVoxelLimits *This ,
                                     const EAxis pAxis);

FQUALIFIER
G4bool GPVoxelLimits_IsXLimited(GPVoxelLimits *This );

FQUALIFIER
G4bool GPVoxelLimits_IsYLimited(GPVoxelLimits *This );

FQUALIFIER
G4bool GPVoxelLimits_IsZLimited(GPVoxelLimits *This );

FQUALIFIER
G4bool GPVoxelLimits_IsLimited(GPVoxelLimits *This );

FQUALIFIER
G4bool GPVoxelLimits_IsLimited2(GPVoxelLimits *This ,
                                const EAxis pAxis);

FQUALIFIER
G4bool GPVoxelLimits_Inside( GPVoxelLimits *This , 
                             GPThreeVector pVec);

}

#endif
