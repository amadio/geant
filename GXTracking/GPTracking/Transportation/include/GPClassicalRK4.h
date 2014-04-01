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
// $Id: G4MagIntegratorStepper.hh,v 1.14 2009-11-05 18:31:15 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4MagIntegratorStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field
//
//  A Stepper must integrate over                NumberOfVariables elements,
//   and also copy (from input to output) any of NoStateVariables  
//   not included in the NumberOfVariables.  
// 
//  So it is expected that NoStateVariables >= NumberOfVariables

// History:
// - 15.01.97  J. Apostolakis (J.Apostolakis@cern.ch)
// --------------------------------------------------------------------

#ifndef GPClassicalRK4_HH
#define GPClassicalRK4_HH

#include "GPTypeDef.h"
#include "GPThreeVector.h"
#include "GPEquationOfMotion.h"
#include "GPClassicalRK4.h"

struct GPClassicalRK4
{
  //-------------------------------------------------
  // class G4ClassicalRK4 : public G4MagErrorStepper 
  //-------------------------------------------------  

  //  G4double *dydxm, *dydxt, *yt; // scratch space - not state 
  G4double dydxm[12], dydxt[12], yt[12]; // scratch space - not state 

  //--------------------------------------------------------
  // class G4MagErrorStepper : public G4MagIntegratorStepper
  //--------------------------------------------------------  

  // STATE
  GPThreeVector fInitialPoint, fMidPoint, fFinalPoint;
  // Data stored in order to find the chord
  
  // Dependent Objects, owned --- part of the STATE 
  //  G4double *yInitial, *yMiddle, *dydxMid, *yOneStep;
  G4double yInitial[12], yMiddle[12], dydxMid[12], yOneStep[12];
  // The following arrays are used only for temporary storage
  // they are allocated at the class level only for efficiency -
  // so that calls to new and delete are not made in Stepper().

  //--------------------------------------------------------
  // class G4MagIntegratorStepper
  //--------------------------------------------------------  

  GPEquationOfMotion *fEquation_Rhs;
  G4int  fNoIntegrationVariables; // Number of Variables in integration
  G4int  fNoStateVariables;       // Number required for FieldTrack

};

extern "C" {

FQUALIFIER
void GPClassicalRK4_Constructor( GPClassicalRK4 *This, 
				 GPEquationOfMotion* EqRhs,
				 G4int numberOfVariables);

FQUALIFIER
void GPClassicalRK4_Destructor( GPClassicalRK4 *This);

FQUALIFIER
G4int GPClassicalRK4_IntegratorOrder();

FQUALIFIER
void GPClassicalRK4_DumbStepper( GPClassicalRK4 *This,
				 const G4double  yIn[],
				 const G4double  dydx[],
				 G4double  h,
				 G4double  yOut[]);

FQUALIFIER
void GPClassicalRK4_G4MagErrorStepper_Constructor(GPClassicalRK4 *This,
						  GPEquationOfMotion *EquationRhs,
						  G4int numberOfVariables, 
						  G4int numStateVariables); 
  
FQUALIFIER
void GPClassicalRK4_G4MagErrorStepper_Destructor(GPClassicalRK4 *This);
  
FQUALIFIER
void GPClassicalRK4_Stepper( GPClassicalRK4 *This, 
			     const G4double yInput[],
			     const G4double dydx[],
			     G4double hstep,
			     G4double yOutput[],
			     G4double yError [] );

FQUALIFIER
  G4double GPClassicalRK4_DistChord(GPClassicalRK4 *This);

FQUALIFIER
void GPClassicalRK4_G4MagIntegratorStepper_Constructor(GPClassicalRK4 *This,
						       GPEquationOfMotion* Equation,
						       G4int       num_integration_vars,
						       G4int       num_state_vars);
  
FQUALIFIER
void GPClassicalRK4_ComputeRightHandSide( GPClassicalRK4 *This,
					  const G4double y[], 
					  G4double dydx[] ) ;

FQUALIFIER
GPEquationOfMotion* GPClassicalRK4_GetEquationOfMotion(GPClassicalRK4 *This);

FQUALIFIER
void GPClassicalRK4_SetEquationOfMotion(GPClassicalRK4 *This,
					GPEquationOfMotion* newEquation);

FQUALIFIER
G4int GPClassicalRK4_GetNumberOfVariables(GPClassicalRK4 *This); 

FQUALIFIER
G4int GPClassicalRK4_GetNumberOfStateVariables(GPClassicalRK4 *This); 

FQUALIFIER 
void GPClassicalRK4_RightHandSide(GPClassicalRK4 *This,
				  const  G4double y[], 
				  G4double dydx[] );   
FQUALIFIER 
void GPClassicalRK4_NormaliseTangentVector(GPClassicalRK4 *This, 
					   G4double vec[6] );

FQUALIFIER 
void GPClassicalRK4_NormalisePolarizationVector(GPClassicalRK4 *This, 
						G4double vec[12] );

}

#endif
