#ifndef GPLOGICALVOLUME_HH
#define GPLOGICALVOLUME_HH

#include "GPGeomdefs.h"
#include "GPMaterial.h"
#include "GPVSolid.h"
#include "GPVPhysicalVolume.h"
#include "GPSmartVoxelHeader.h"

struct GPLogicalVolume
{
  // unique index, can be used for references and for sharing with
  // another representation of the geometry.
  // -1 indicates it has not yet been set.
  int fIndex;
  G4int fNoDaughters;
  struct GPVPhysicalVolume* GEOMETRYLOC *fDaughters;
  GEOMETRYLOC GPMaterial* fMaterial;
  GEOMETRYLOC struct GPVSolid* fSolid;
  GEOMETRYLOC struct GPSmartVoxelHeader *fVoxel;
};

extern "C" {

FQUALIFIER
void GPLogicalVolume_Constructor(GPLogicalVolume *This, 
				 GPVSolid* pSolid, 
				 GPMaterial* pMaterial );

FQUALIFIER
G4int GPLogicalVolume_GetNoDaughters(const GPLogicalVolume* This);

FQUALIFIER
GPVPhysicalVolume* GPLogicalVolume_GetDaughter(const GPLogicalVolume* This, 
					       const G4int i);

FQUALIFIER
void GPLogicalVolume_AddDaughter(GPLogicalVolume* This, 
				 GPVPhysicalVolume* pNewDaughter);

FQUALIFIER
GPVSolid* GPLogicalVolume_GetSolid(const GPLogicalVolume* This);

FQUALIFIER
void GPLogicalVolume_SetSolid(GPLogicalVolume *This,
			      GPVSolid *pSolid);

FQUALIFIER
GPMaterial* GPLogicalVolume_GetMaterial(const GPLogicalVolume* This);

FQUALIFIER
void GPLogicalVolume_SetMaterial(GPLogicalVolume *This,
				 GPMaterial *pMaterial);  

FQUALIFIER
void GPLogicalVolume_UpdateMaterial(GPLogicalVolume *This,
				    GPMaterial *pMaterial);

FQUALIFIER
GPSmartVoxelHeader* GPLogicalVolume_GetVoxelHeader(GPLogicalVolume* This);

FQUALIFIER
void GPLogicalVolume_SetVoxelHeader(GPLogicalVolume* This, 
				    GPSmartVoxelHeader * pVoxel);
}

#endif
