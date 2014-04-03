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
// $Id: G4CashKarpRKF45.hh,v 1.12 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4CashKarpRKF45
//
// Class description:
//
// The Cash-Karp Runge-Kutta-Fehlberg 4/5 method is an embedded fourth
// order method (giving fifth-order accuracy) for the solution of an ODE.
// Two different fourth order estimates are calculated; their difference
// gives an error estimate. [ref. Numerical Recipes in C, 2nd Edition]
// It is used to integrate the equations of the motion of a particle 
// in a magnetic field.

// History:
// - Created. J.Apostolakis, V.Grichine - 30.1.97
// -------------------------------------------------------------------

#ifndef GPCashKarpRKF45_HH
#define GPCashKarpRKF45_HH

#include "GPTypeDef.h"
#include "GPEquationOfMotion.h"
#include "GPCashKarpRKF45.h"
#include "GPLineSection.h"

struct GPCashKarpRKF45
{
  //--------------------------------------------------------
  // class GPCashKarpRKF45 : public G4MagIntegratorStepper
  //--------------------------------------------------------  

  //  G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *yTemp, *yIn;
  // scratch space
  G4double ak2[12], ak3[12], ak4[12], ak5[12], ak6[12], ak7[12];
  G4double yTemp[12],yIn[12];

  G4double fLastStepLength;
  G4double fLastInitialVector[12], fLastFinalVector[12],
    fLastDyDx[12], fMidVector[12], fMidError[12];
  
  //--------------------------------------------------------
  // class G4MagIntegratorStepper
  //--------------------------------------------------------  

  GPEquationOfMotion *fEquation_Rhs;
  G4int  fNoIntegrationVariables; // Number of Variables in integration
  G4int  fNoStateVariables;       // Number required for FieldTrack

};

extern "C" {

FQUALIFIER
void GPCashKarpRKF45_Constructor( GPCashKarpRKF45 *This, 
				  GPEquationOfMotion* EqRhs, 
				  G4int numberOfVariables,
				  G4bool primary);

FQUALIFIER
G4int GPCashKarpRKF45_IntegratorOrder();

FQUALIFIER
void GPCashKarpRKF45_Stepper( GPCashKarpRKF45 *This, 
			      const G4double yInput[],
			      const G4double dydx[],
			      G4double Step,
			      G4double yOut[],
			      G4double yErr[]);

FQUALIFIER
G4double GPCashKarpRKF45_DistChord(GPCashKarpRKF45 *This);

FQUALIFIER
void GPCashKarpRKF45_G4MagIntegratorStepper_Constructor(GPCashKarpRKF45 *This,
							GPEquationOfMotion* Equation,
							G4int       num_integration_vars,
							G4int       num_state_vars);

FQUALIFIER
void GPCashKarpRKF45_ComputeRightHandSide( GPCashKarpRKF45 *This,
					   const G4double y[],
					   G4double dydx[] ) ;

FQUALIFIER
GPEquationOfMotion* GPCashKarpRKF45_GetEquationOfMotion(GPCashKarpRKF45 *This);

FQUALIFIER
void GPCashKarpRKF45_SetEquationOfMotion(GPCashKarpRKF45 *This,
					 GPEquationOfMotion* newEquation);

FQUALIFIER
G4int GPCashKarpRKF45_GetNumberOfVariables(GPCashKarpRKF45 *This);

FQUALIFIER
G4int GPCashKarpRKF45_GetNumberOfStateVariables(GPCashKarpRKF45 *This);

FQUALIFIER 
void GPCashKarpRKF45_RightHandSide(GPCashKarpRKF45 *This,
				   const  G4double y[], 
				   G4double dydx[] );

FQUALIFIER 
void GPCashKarpRKF45_NormaliseTangentVector(GPCashKarpRKF45 *This, 
					    G4double vec[6] );

FQUALIFIER 
void GPCashKarpRKF45_NormalisePolarizationVector(GPCashKarpRKF45 *This, 
						 G4double vec[12] );


}

#endif
