#include "GPVPhysicalVolume.h"

FQUALIFIER
void GPVPhysicalVolume_Constructor(GPVPhysicalVolume *This,
				   GPRotationMatrix pRot,
				   GPThreeVector tlate,
				   GPLogicalVolume* pLogical )
{
  This->fIndex = -1;
  This->frot = pRot;
  This->ftrans = tlate;
  This->flogical = pLogical;
  This->flmother = NULL;
}

FQUALIFIER
G4int GPVPhysicalVolume_GetMultiplicity(GEOMETRYLOC GPVPhysicalVolume *This )
{
  return 1;
}

FQUALIFIER
GPRotationMatrix* GPVPhysicalVolume_GetObjectRotation(GEOMETRYLOC GPVPhysicalVolume *This )
{
  GPRotationMatrix  aRotM; 
  GPRotationMatrix  IdentityRM = GPRotationMatrix_create(1,0,0,
							 0,1,0,
							 0,0,1);
  GPRotationMatrix* retval; 

  // Insure against frot being a null pointer
  if(&(This->frot))
  {
    aRotM = GPRotationMatrix_inverse(&(This->frot));
    retval= &aRotM;
  }
  else
  {
    retval= &IdentityRM;
  }
  return retval;
}

FQUALIFIER
G4bool GPVPhysicalVolume_CheckOverlaps(G4int, G4double, G4bool)
{
  return false;
}

FQUALIFIER
GPThreeVector GPVPhysicalVolume_GetTranslation(GEOMETRYLOC GPVPhysicalVolume *This)
{
  return This->ftrans;
}

FQUALIFIER
void GPVPhysicalVolume_SetTranslation(GPVPhysicalVolume *This,
				      GPThreeVector v)
{
  This->ftrans = v;
}

FQUALIFIER
GPRotationMatrix* GPVPhysicalVolume_GetRotation(GEOMETRYLOC GPVPhysicalVolume *This )
{
  return &(This->frot);
}

FQUALIFIER
void GPVPhysicalVolume_SetRotation( GPVPhysicalVolume *This,
				    GPRotationMatrix pRot)
{
  This->frot=pRot;
}

FQUALIFIER
GEOMETRYLOC GPLogicalVolume* GPVPhysicalVolume_GetLogicalVolume(GEOMETRYLOC GPVPhysicalVolume *This)
{
  return This->flogical;
}

FQUALIFIER
void GPVPhysicalVolume_SetLogicalVolume(GPVPhysicalVolume *This,
					GPLogicalVolume *pLogical)
{
  This->flogical=pLogical;
}

FQUALIFIER
GEOMETRYLOC GPLogicalVolume* GPVPhysicalVolume_GetMotherLogical(GEOMETRYLOC GPVPhysicalVolume *This)
{
  return This->flmother;
}

FQUALIFIER
void GPVPhysicalVolume_SetMotherLogical(GPVPhysicalVolume *This,
					GPLogicalVolume *pMother)
{
  This->flmother=pMother;
}

FQUALIFIER
GPRotationMatrix GPVPhysicalVolume_GetObjectRotationValue(GEOMETRYLOC GPVPhysicalVolume *This)
{
  GPRotationMatrix  aRotM;   // Initialised to identity

  // Insure against frot being a null pointer
  if(&(This->frot))
  {
    aRotM =  GPRotationMatrix_inverse(&(This->frot));
  }
  return aRotM;
}

FQUALIFIER
GPThreeVector  GPVPhysicalVolume_GetObjectTranslation(GEOMETRYLOC GPVPhysicalVolume *This)
{
  return  This->ftrans;
}

FQUALIFIER
GPRotationMatrix* GPVPhysicalVolume_GetFrameRotation(GEOMETRYLOC GPVPhysicalVolume *This)
{
  return &(This->frot);
}

FQUALIFIER
GPThreeVector GPVPhysicalVolume_GetFrameTranslation(GEOMETRYLOC GPVPhysicalVolume *This)
{
  return GPThreeVector_mult(This->ftrans,-1.0);
}
