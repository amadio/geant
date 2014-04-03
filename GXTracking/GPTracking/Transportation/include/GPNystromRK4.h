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
// $Id: G4NystromRK4.hh,v 1.4 2010-07-14 10:00:36 gcosmo Exp $ 
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4NystromRK4
//
// Class description:
//
// Integrate the equations of the motion of a particle in a magnetic field
// using 4th Runge-Kutta-Nystrom method with errors estimation 
// (ATL-SOFT-PUB-2009-01)
// Current form can be used only for 'pure' magnetic field.
// Notes: 1) field must be time-independent.
//        2) time is not integrated
// 
// History:
// - Created: I.Gavrilenko   15.05.2009   (as G4AtlasRK4)
// - Adaptations:  J. Apostolakis  May-Nov 2009
// -------------------------------------------------------------------

#ifndef GPNYSTROMRK4_HH
#define GPNYSTROMRK4_HH

//#include "globals.hh"
//#include "G4MagIntegratorStepper.hh"
//#include "G4Mag_EqRhs.hh"

#include "GPTypeDef.h"
#include "GPEquationOfMotion.h"

struct GPNystromRK4
{
  //--------------------------------------------------------
  // class GPNystromRK4 : public G4MagIntegratorStepper
  //--------------------------------------------------------

    G4double      m_lastField[3];
    G4double      m_fldPosition[4];
    G4double      m_magdistance ;
    G4double      m_magdistance2;
    G4double      m_cof         ;
    G4double      m_mom         ;
    G4double      m_imom        ;
    G4bool        m_cachedMom   ;
    G4double      m_iPoint   [3];
    G4double      m_mPoint   [3];
    G4double      m_fPoint   [3];

  //--------------------------------------------------------
  // class G4MagIntegratorStepper
  //--------------------------------------------------------  

  GPEquationOfMotion *m_fEq;
  G4int  fNoIntegrationVariables; // Number of Variables in integration
  G4int  fNoStateVariables;       // Number required for FieldTrack
    
};

extern "C" {

FQUALIFIER
void GPNystromRK4_Constructor(GPNystromRK4 *This, 
			      GPEquationOfMotion* EqRhs, 
			      G4double distanceConstField);

FQUALIFIER
void GPNystromRK4_Stepper(GPNystromRK4 *This,
			  const G4double P[],const G4double dPdS[],
			  G4double Step,G4double Po[],G4double Err[]);

FQUALIFIER
G4double GPNystromRK4_DistChord(GPNystromRK4 *This);

FQUALIFIER
void  GPNystromRK4_ComputeRightHandSide(GPNystromRK4 *This, 
					const G4double P[], 
					G4double dPdS[]);

FQUALIFIER
void  GPNystromRK4_SetDistanceForConstantField(GPNystromRK4 *This, 
					       G4double length );

FQUALIFIER
G4double  GPNystromRK4_GetDistanceForConstantField(GPNystromRK4 *This);

FQUALIFIER
void GPNystromRK4_getField (GPNystromRK4 *This, const G4double P[4]);

FQUALIFIER
void GPNystromRK4_G4MagIntegratorStepper_Constructor(GPNystromRK4 *This,
                                         GPEquationOfMotion* Equation,
                                         G4int       num_integration_vars,
                                         G4int       num_state_vars);
}

#endif  // GPNYSTROMRK4
