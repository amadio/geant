#ifndef GPVPHYSICALVOLUME_HH
#define GPVPHYSICALVOLUME_HH

#include "GPGeomdefs.h"
#include "GPRotationMatrix.h"
#include "GPLogicalVolume.h"

struct GPVPhysicalVolume
{
   // unique index, can be used for references and for sharing with
   // another representation of the geometry.
   // -1 indicates it has not yet been set.
   int fIndex;
   GPRotationMatrix frot;
   GPThreeVector ftrans;
   GEOMETRYLOC struct GPLogicalVolume *flogical;
   GEOMETRYLOC struct GPLogicalVolume *flmother;
};

extern "C" {

  FQUALIFIER
  void GPVPhysicalVolume_Constructor(
                GPVPhysicalVolume *This,
                GPRotationMatrix pRot,
                GPThreeVector tlate,
                GPLogicalVolume* pLogical );


FQUALIFIER
G4int GPVPhysicalVolume_GetMultiplicity(GEOMETRYLOC GPVPhysicalVolume *This );

FQUALIFIER
GPRotationMatrix* GPVPhysicalVolume_GetObjectRotation(GEOMETRYLOC GPVPhysicalVolume *This);

FQUALIFIER
G4bool GPVPhysicalVolume_CheckOverlaps(G4int, G4double, G4bool);

FQUALIFIER
GPThreeVector GPVPhysicalVolume_GetTranslation(GEOMETRYLOC GPVPhysicalVolume *This);

FQUALIFIER
void GPVPhysicalVolume_SetTranslation(GPVPhysicalVolume *This,
				      GPThreeVector v);

FQUALIFIER
GPRotationMatrix* GPVPhysicalVolume_GetRotation(GEOMETRYLOC GPVPhysicalVolume *This );

FQUALIFIER
void GPVPhysicalVolume_SetRotation( GPVPhysicalVolume *This,
				    GPRotationMatrix pRot);

FQUALIFIER
GEOMETRYLOC GPLogicalVolume* GPVPhysicalVolume_GetLogicalVolume(GEOMETRYLOC GPVPhysicalVolume *This);

FQUALIFIER
void GPVPhysicalVolume_SetLogicalVolume(GPVPhysicalVolume *This,
					GPLogicalVolume *pLogical);

FQUALIFIER
GEOMETRYLOC GPLogicalVolume* GPVPhysicalVolume_GetMotherLogical(GEOMETRYLOC GPVPhysicalVolume *This);

  FQUALIFIER
  void GPVPhysicalVolume_SetMotherLogical(GPVPhysicalVolume *This, 
					  GPLogicalVolume *pMother);

FQUALIFIER
GPRotationMatrix GPVPhysicalVolume_GetObjectRotationValue(GEOMETRYLOC GPVPhysicalVolume *This);

FQUALIFIER
GPThreeVector  GPVPhysicalVolume_GetObjectTranslation(GEOMETRYLOC GPVPhysicalVolume *This);

FQUALIFIER
GPRotationMatrix* GPVPhysicalVolume_GetFrameRotation(GEOMETRYLOC GPVPhysicalVolume *This);

FQUALIFIER
GPThreeVector GPVPhysicalVolume_GetFrameTranslation(GEOMETRYLOC GPVPhysicalVolume *This);

}

#endif
