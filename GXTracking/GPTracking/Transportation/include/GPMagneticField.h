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
// $Id: G4MagneticField.hh,v 1.14 2006-06-29 18:23:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4MagneticField
//
// Class description:
//
// Magnetic Field abstract class, implements inquiry function interface.

// History:
// - Created. JA, January 13th, 1996.
// --------------------------------------------------------------------

#ifndef GPMAGNETICFIELD_HH
#define GPMAGNETICFIELD_HH

#include "GPTypeDef.h"
#include "GPFieldMap.h"
#include "GPMagneticField.h"

struct GPMagneticField
{
  G4bool  fGravityActive;
  GPFieldMap *fieldMap;

  G4double prm[9];
  G4double hb0;
  G4double hlova;
  G4double ainv;
  G4double coeff;

};

extern "C" {

FQUALIFIER 
void GPMagneticField_Constructor(GPMagneticField *This,
                                 GPFieldMap *fieldMapArray);

FQUALIFIER 
void GPMagneticField_Constructor_Parametric(GPMagneticField *This);

FQUALIFIER 
void  GPMagneticField_GetFieldValue( GPMagneticField *This,
                                     const  G4double point[4],
                                     double *bField );

FQUALIFIER 
void  GPMagneticField_GetFieldValue_Parametric( GPMagneticField *This,
                                                const  G4double point[4],
                                                G4double *bField );

FQUALIFIER
void GPMagneticField_ffunkti(G4double u, G4double *ff);

FQUALIFIER
void  GPMagneticField_SetGravityActive(GPMagneticField *This, 
                                       G4bool OnOffFlag );

FQUALIFIER
G4bool GPMagneticField_IsGravityActive(GPMagneticField *This);

FQUALIFIER
G4bool GPMagneticField_DoesFieldChangeEnergy(GPMagneticField *This);

}

#endif
