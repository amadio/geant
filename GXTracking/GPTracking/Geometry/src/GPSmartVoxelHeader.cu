
// $Id: GPSmartVoxelHeader.icc,v 1.6 2006-06-29 18:32:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GPSmartVoxelHeader Inline implementation
//
// --------------------------------------------------------------------

#include "GPSmartVoxelHeader.h"

FQUALIFIER
G4int GPSmartVoxelHeader_GetMaxEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelHeader *This)
{
  return This->fmaxEquivalent;
}

FQUALIFIER
void GPSmartVoxelHeader_SetMaxEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelHeader *This,
						G4int pMax)
{
  This->fmaxEquivalent=pMax;
}

FQUALIFIER
G4int GPSmartVoxelHeader_GetMinEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelHeader *This)
{
  return This->fminEquivalent;
}

FQUALIFIER
void GPSmartVoxelHeader_SetMinEquivalentSliceNo(GEOMETRYLOC GPSmartVoxelHeader *This,
						G4int pMin)
{
  This->fminEquivalent=pMin;
}

FQUALIFIER
EAxis GPSmartVoxelHeader_GetAxis(GEOMETRYLOC GPSmartVoxelHeader *This)
{
  return This->faxis;
}

FQUALIFIER
EAxis GPSmartVoxelHeader_GetParamAxis(GEOMETRYLOC GPSmartVoxelHeader *This)
{
  return This->fparamAxis;
}

FQUALIFIER
G4double GPSmartVoxelHeader_GetMaxExtent(GEOMETRYLOC GPSmartVoxelHeader *This)
{
  return This->fmaxExtent;
}

FQUALIFIER
G4double GPSmartVoxelHeader_GetMinExtent(GEOMETRYLOC GPSmartVoxelHeader *This)
{
  return This->fminExtent;
}

FQUALIFIER
G4int GPSmartVoxelHeader_GetNoSlices(GEOMETRYLOC GPSmartVoxelHeader *This)
{
  return This->fnumSlices;
}

FQUALIFIER
GPSmartVoxelProxy* GPSmartVoxelHeader_GetSlice(GEOMETRYLOC GPSmartVoxelHeader *This,
					       G4int n)
{
  return This->fslices[n];
}
