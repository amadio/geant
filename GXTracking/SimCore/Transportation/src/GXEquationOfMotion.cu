#include "GXEquationOfMotion.h"
#include "GXMagneticField.h"
#include "GPThreeVector.h"

FQUALIFIER
void GXEquationOfMotion_Constructor( GXEquationOfMotion *This,
			             GXMagneticField* magField,
                                     G4double charge)
{
  This->fField = magField;
  GXEquationOfMotion_SetCharge(This,charge);
}


FQUALIFIER
void GXEquationOfMotion_RightHandSide( GXEquationOfMotion *This, 
				       const  G4double y[],
				       G4double dydx[]   )
{
  // Position
  G4double Position[3] = {y[0], y[1], y[2]};
  
  // B-field
  G4double B[3];   
  GXMagneticField_GetFieldValue(This->fField, Position, B);

  GPThreeVector p = {y[3], y[4], y[5]};
  G4double momentum_mag_square = p.x*p.x + p.y*p.y + p.z*p.z;
  G4double inv_momentum_magnitude = 1.0 / sqrt( momentum_mag_square );
  
  G4double cof = This->fCoeff*inv_momentum_magnitude;

  //the right hand side of the equation of motion
  
  dydx[0] = p.x*inv_momentum_magnitude;    // (d/ds)x = Vx/V
  dydx[1] = p.y*inv_momentum_magnitude;    // (d/ds)y = Vy/V
  dydx[2] = p.z*inv_momentum_magnitude;    // (d/ds)z = Vz/V
  
  dydx[3] = cof*(p.y*B[2] - p.z*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
  dydx[4] = cof*(p.z*B[0] - p.x*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
  dydx[5] = cof*(p.x*B[1] - p.y*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)
}


FQUALIFIER
void GXEquationOfMotion_SetCharge( GXEquationOfMotion *This,
                                   G4double charge)
{
  This->fCharge = charge; 
  This->fCoeff = charge*299.792458; //charge in eplus unit times c (c_light) 
}

FQUALIFIER
G4double GXEquationOfMotion_GetCoeff( GXEquationOfMotion *This)
{
  return This->fCoeff;
}

FQUALIFIER
void GXEquationOfMotion_GetFieldValue( GXEquationOfMotion *This,
                                       const  G4double Point[3],
                                       G4double Field[] )
{
  GXMagneticField_GetFieldValue(This->fField, Point, Field );
}
