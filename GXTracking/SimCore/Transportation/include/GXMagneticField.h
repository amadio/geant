#ifndef GXMAGNETICFIELD_HH
#define GXMAGNETICFIELD_HH

#include "GPTypeDef.h"
#include "GXFieldMap.h"

struct GXMagneticField
{
  GXFieldMap *fieldMap;
};

extern "C" {

FQUALIFIER 
void GXMagneticField_Constructor( GXMagneticField *This,
                                  GXFieldMap *fieldMapArray);

FQUALIFIER 
void GXMagneticField_GetFieldValue( GXMagneticField *This,
				    const  G4double point[3],
				    G4double *bField );
}

#endif
