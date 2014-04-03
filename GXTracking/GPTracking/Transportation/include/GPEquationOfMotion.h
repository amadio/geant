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
// $Id: G4EquationOfMotion.hh,v 1.10 2006-06-29 18:22:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4EquationOfMotion
//
// Class description:
//
// Abstract Base Class for the right hand size of the equation of
// motion of a particle in a field.

// History:
// - Created. J.Apostolakis
// -------------------------------------------------------------------

#ifndef GPEquationOfMotion_HH
#define GPEquationOfMotion_HH

#include "GPTypeDef.h"        // "globals.hh"
#include "GPMagneticField.h"  // required in inline method implementations
#include "GPEquationOfMotion.h"

//this is the based class of G4Mag_EqRhs, G4Mag_UsualEqRhs

// static const int G4maximum_number_of_field_components = 16;
enum { G4maximum_number_of_field_components = 16 } ;

struct GPEquationOfMotion 
{
  GPMagneticField *itsField;
  G4double fCof_val;  //G4Mag_EqRhs
  G4double fInvCurrentMomentumXc;
  G4double fUnitConstant; //G4Mag_EqRhs
};

extern "C" {

FQUALIFIER
void GPEquationOfMotion_Constructor( GPEquationOfMotion *This, 
			             GPMagneticField* magField);

FQUALIFIER
GPMagneticField* GPEquationOfMotion_GetFieldObj( GPEquationOfMotion *This );

FQUALIFIER
void GPEquationOfMotion_SetFieldObj( GPEquationOfMotion *This,
				     GPMagneticField* pField);

FQUALIFIER
void GPEquationOfMotion_GetFieldValue( GPEquationOfMotion *This,
				       const  G4double Point[4],
				       G4double Field[] ) ;

FQUALIFIER
void GPEquationOfMotion_RightHandSide( GPEquationOfMotion *This,
				       const  G4double y[],
				       G4double dydx[]   );

FQUALIFIER
void GPEquationOfMotion_EvaluateRhsReturnB(GPEquationOfMotion *This,
					   const G4double y[],
					   G4double dydx[],
					   G4double  Field[]  );

FQUALIFIER
void GPEquationOfMotion_EvaluateRhsGivenB( GPEquationOfMotion *This,
					   const G4double y[],
					   const G4double B[3],
                                           G4double dydx[] );

FQUALIFIER
G4double GPEquationOfMotion_FCof(GPEquationOfMotion *This);

FQUALIFIER
void GPEquationOfMotion_SetChargeMomentumMass( GPEquationOfMotion *This,
					       G4double particleCharge, // e+ units
					       G4double MomentumXc,
					       G4double mass);

}

#endif
