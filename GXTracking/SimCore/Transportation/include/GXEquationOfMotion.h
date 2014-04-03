#ifndef GXEquationOfMotion_HH
#define GXEquationOfMotion_HH

#include "GPTypeDef.h"       
#include "GXMagneticField.h"

struct GXEquationOfMotion 
{
  GXMagneticField *fField;
  G4double fCharge; 
  G4double fCoeff; 
};

extern "C" {

FQUALIFIER
void GXEquationOfMotion_Constructor( GXEquationOfMotion *This, 
			             GXMagneticField* magField,
				     G4double charge);

FQUALIFIER
void GXEquationOfMotion_RightHandSide( GXEquationOfMotion *This,
				       const  G4double y[],
				       G4double dydx[]   );

FQUALIFIER
void GXEquationOfMotion_SetCharge( GXEquationOfMotion *This,
				   G4double charge);

FQUALIFIER
G4double GXEquationOfMotion_GetCoeff( GXEquationOfMotion *This);

FQUALIFIER
void GXEquationOfMotion_GetFieldValue( GXEquationOfMotion *This,
                                       const  G4double Point[3],
                                       G4double Field[] );
}

#endif
