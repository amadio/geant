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
// $Id: G4EquationOfMotion.cc,v 1.9 2006-06-29 18:23:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

#include "GPEquationOfMotion.h"
#include "GPMagneticField.h"
#include "GPConstants.h"

FQUALIFIER
void GPEquationOfMotion_Constructor( GPEquationOfMotion *This, 
			             GPMagneticField* magField)
{
  This->itsField = magField;
  This->fCof_val = 0.0;
  This->fInvCurrentMomentumXc = 0.0;
  This->fUnitConstant = 0.299792458 * (GeV/tesla*meter); //(GeV/(tesla*m))
}

FQUALIFIER
GPMagneticField* GPEquationOfMotion_GetFieldObj( GPEquationOfMotion *This )
{
  return This->itsField;
}

FQUALIFIER
void GPEquationOfMotion_SetFieldObj( GPEquationOfMotion *This,
				     GPMagneticField* pField)
{
  This->itsField= pField;
}

FQUALIFIER
void GPEquationOfMotion_GetFieldValue( GPEquationOfMotion *This,
				       const  G4double Point[4],
				       G4double Field[] ) 
{
  GPMagneticField_GetFieldValue(This->itsField, Point, Field );
}

FQUALIFIER
void GPEquationOfMotion_RightHandSide( GPEquationOfMotion *This,
				  const  G4double y[],
				  G4double dydx[]   )
{
  //  G4double Field[G4maximum_number_of_field_components];   
  G4double B[4];   
  G4double  PositionAndTime[4];
  
  // Position
  PositionAndTime[0] = y[0];
  PositionAndTime[1] = y[1];
  PositionAndTime[2] = y[2];
  // Global Time
  PositionAndTime[3] = y[7];  // See G4FieldTrack::LoadFromArray
  
  //  GPEquationOfMotion_GetFieldValue(This,PositionAndTime, Field) ;
  GPEquationOfMotion_GetFieldValue(This,PositionAndTime, B) ;

  //  GPEquationOfMotion_EvaluateRhsGivenB(This, y, Field, dydx );

  G4double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
  G4double inv_momentum_magnitude = 1.0 / sqrt( momentum_mag_square );
  
  G4double cof = GPEquationOfMotion_FCof(This)*inv_momentum_magnitude;
  
  dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
  dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
  dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V
  
  dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
  dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
  dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)
  
}

FQUALIFIER
void GPEquationOfMotion_EvaluateRhsGivenB( GPEquationOfMotion *This,
					   const G4double y[],
					   const G4double B[3],
                                           G4double dydx[] )
{
   G4double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
   G4double inv_momentum_magnitude = 1.0 / sqrt( momentum_mag_square );

   G4double cof = GPEquationOfMotion_FCof(This)*inv_momentum_magnitude;

   dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
   dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
   dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

   dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)

}

FQUALIFIER
G4double GPEquationOfMotion_FCof(GPEquationOfMotion *This)
{
  return This->fCof_val;
}

FQUALIFIER
void GPEquationOfMotion_SetChargeMomentumMass( GPEquationOfMotion *This,
					       G4double particleCharge, // e+ units
					       G4double MomentumXc,
					       G4double mass)
{
   This->fInvCurrentMomentumXc =1.0 / MomentumXc;
   This->fCof_val = particleCharge*eplus*c_light ; //  B must be in Tesla
   //fCof_val = fUnitConstant*particleCharge/MomentumXc; //  B must be in Tesla
   //fMass = particleMass;
}


FQUALIFIER
void GPEquationOfMotion_EvaluateRhsReturnB(GPEquationOfMotion *This,
					   const G4double y[],
					   G4double dydx[],
					   G4double  Field[]  )
{
  G4double  PositionAndTime[4];
  
  // Position
  PositionAndTime[0] = y[0];
  PositionAndTime[1] = y[1];
  PositionAndTime[2] = y[2];
  // Global Time
  PositionAndTime[3] = y[7];  // See G4FieldTrack::LoadFromArray
  
  GPEquationOfMotion_GetFieldValue(This,PositionAndTime, Field) ;
  GPEquationOfMotion_EvaluateRhsGivenB(This, y, Field, dydx );
}


