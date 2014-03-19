#include "GPLogicalVolume.h"

FQUALIFIER
void GPLogicalVolume_Constructor(GPLogicalVolume *This, 
				 GPVSolid* pSolid, 
				 GPMaterial* pMaterial )
{
  This->fIndex = 0;
  This->fNoDaughters = 0;
  This->fDaughters = NULL;
  This->fMaterial = pMaterial;
  This->fSolid = pSolid;
  This->fVoxel = 0;
}

FQUALIFIER
G4int GPLogicalVolume_GetNoDaughters(GEOMETRYLOC const GPLogicalVolume* This)
{
  return This->fNoDaughters;
}

FQUALIFIER
GEOMETRYLOC GPVPhysicalVolume* GPLogicalVolume_GetDaughter(GEOMETRYLOC const GPLogicalVolume* This, 
							   const G4int i)
{
  return This->fDaughters[i];
}

FQUALIFIER
void GPLogicalVolume_AddDaughter(GPLogicalVolume* This, 
				 GPVPhysicalVolume* pNewDaughter)
{
  This->fNoDaughters++;
  This->fDaughters[This->fNoDaughters-1] = pNewDaughter;
}

FQUALIFIER
GEOMETRYLOC GPVSolid* GPLogicalVolume_GetSolid(GEOMETRYLOC const GPLogicalVolume* This)
{
  return This->fSolid;
}

FQUALIFIER
void GPLogicalVolume_SetSolid(GPLogicalVolume *This,
			      GPVSolid *pSolid)
{
  This->fSolid=pSolid;
}

FQUALIFIER
GEOMETRYLOC GPMaterial* GPLogicalVolume_GetMaterial(GEOMETRYLOC const GPLogicalVolume* This)
{
  return This->fMaterial;
}

FQUALIFIER
void GPLogicalVolume_SetMaterial(GPLogicalVolume *This,
				 GPMaterial *pMaterial)
{
  This->fMaterial=pMaterial;
}

FQUALIFIER
void GPLogicalVolume_UpdateMaterial(GPLogicalVolume *This,
				    GPMaterial *pMaterial)
{
  This->fMaterial=pMaterial;
}

GEOMETRYLOC GPSmartVoxelHeader* GPLogicalVolume_GetVoxelHeader(GEOMETRYLOC GPLogicalVolume* This)
{
  return This->fVoxel;
}

void GPLogicalVolume_SetVoxelHeader(GPLogicalVolume* This, 
				    GPSmartVoxelHeader * pVoxel)
{
  This->fVoxel = pVoxel;
}
