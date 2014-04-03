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
// $Id: G4MagIntegratorDriver.hh,v 1.21 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4MagInt_Driver
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and insures that 
// the error is within acceptable bounds.

// History:
// - Created. J.Apostolakis.
// --------------------------------------------------------------------

#ifndef GPMagInt_Driver_HH
#define GPMagInt_Driver_HH

#include "GPTypeDef.h"
#include "GPUtils.h"
#include "GPFieldTrack.h"
#include "GPClassicalRK4.h"
#include "GPMagIntegratorDriver.h"

struct GPMagInt_Driver
{
  // ---------------------------------------------------------------
  //  INVARIANTS 
  
  G4double  fMinimumStep;
  // Minimum Step allowed in a Step (in absolute units)
  G4double  fSmallestFraction;      //   Expected range 1e-12 to 5e-15;  
  // Smallest fraction of (existing) curve length - in relative units
  //  below this fraction the current step will be the last 
  
  G4int  fNoIntegrationVariables;  // Number of Variables in integration
  G4int  fMinNoVars;               // Minimum number for FieldTrack
  G4int  fNoVars;                  // Full number of variable
  
  G4int   fMaxNoSteps;
  //     static const 
  G4int  fMaxStepBase;  
  
  G4double safety;
  G4double pshrnk;   //  exponent for shrinking
  G4double pgrow;    //  exponent for growth
  G4double errcon;
  // Parameters used to grow and shrink trial stepsize.
  
  G4int    fStatisticsVerboseLevel;
  
  // ---------------------------------------------------------------
  // DEPENDENT Objects
  //  GPMagIntegratorStepper *pIntStepper;
  GPClassicalRK4 *pIntStepper;
  
  // ---------------------------------------------------------------
  //  STATE
  
  G4int  fNoTotalSteps, fNoBadSteps, fNoSmallSteps, fNoInitialSmallSteps; 
  G4double fDyerr_max, fDyerr_mx2;
  G4double fDyerrPos_smTot, fDyerrPos_lgTot, fDyerrVel_lgTot; 
  G4double fSumH_sm, fSumH_lg; 
  // Step Statistics 
  
  G4int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
  // Could be varied during tracking - to help identify issues
  
};

//  const G4double max_stepping_increase;
//  const G4double max_stepping_decrease;
//  const G4int  fMaxStepBase;  

//const G4double G4MagInt_Driver::max_stepping_increase = 5.0;
//const G4double G4MagInt_Driver::max_stepping_decrease = 0.1;
const G4double max_stepping_increase = 5.0;
const G4double max_stepping_decrease = 0.1;

//  The (default) maximum number of steps is Base
//  divided by the order of Stepper
//
//const G4int  G4MagInt_Driver::fMaxStepBase = 250;  // Was 5000
const G4int  fMaxStepBase = 250;  // Was 5000


extern "C" {

FQUALIFIER
void GPMagInt_Driver_Constructor(GPMagInt_Driver *This, 
                                 G4double         hminimum, 
                                 GPClassicalRK4  *pStepper,
                                 G4int            numComponents,
                                 G4int            statisticsVerbose);
FQUALIFIER
G4bool GPMagInt_Driver_AccurateAdvance( GPMagInt_Driver *This,
                                        GPFieldTrack& y_current,
                                        G4double     hstep,
                                        G4double     eps,
                                        G4double hinitial );
FQUALIFIER
void GPMagInt_Driver_OneGoodStep( GPMagInt_Driver *This,
                                  G4double y[],      
                                  const G4double dydx[],
                                  G4double& x,       
                                  G4double htry,
                                  G4double eps_rel_max,
                                  G4double& hdid,    
                                  G4double& hnext ); 
FQUALIFIER
G4bool  GPMagInt_Driver_QuickAdvance2(GPMagInt_Driver *This,       
				      GPFieldTrack& y_posvel,
				      const G4double     dydx[],  
				      G4double     hstep,    
				      G4double&    dchord_step,
				      G4double&    dyerr_pos_sq,
				      G4double&    dyerr_mom_rel_sq );  
FQUALIFIER
G4bool  GPMagInt_Driver_QuickAdvance(GPMagInt_Driver *This,        
                                     GPFieldTrack& y_posvel,
                                     const G4double     dydx[],  
                                     G4double     hstep, 
                                     G4double&    dchord_step,
                                     G4double&    dyerr );
FQUALIFIER
G4double GPMagInt_Driver_ComputeNewStepSize(GPMagInt_Driver *This, 
					    G4double  errMaxNorm, 
					    G4double  hstepCurrent);
FQUALIFIER
G4double GPMagInt_Driver_ComputeNewStepSize_WithinLimits(GPMagInt_Driver *This, 
							 G4double errMaxNorm,
							 G4double hstepCurrent);
FQUALIFIER
void GPMagInt_Driver_SetSmallestFraction(GPMagInt_Driver *This,
                                         G4double newFraction);
FQUALIFIER
G4double GPMagInt_Driver_GetHmin( GPMagInt_Driver *This );

FQUALIFIER
G4double GPMagInt_Driver_Hmin( GPMagInt_Driver *This );

FQUALIFIER
G4double GPMagInt_Driver_GetSafety( GPMagInt_Driver *This );

FQUALIFIER
G4double GPMagInt_Driver_GetPshrnk( GPMagInt_Driver *This );

FQUALIFIER
G4double GPMagInt_Driver_GetPgrow( GPMagInt_Driver *This );

FQUALIFIER
G4double GPMagInt_Driver_GetErrcon( GPMagInt_Driver *This );

FQUALIFIER
void GPMagInt_Driver_SetHmin( GPMagInt_Driver *This, G4double newval);

FQUALIFIER
G4double GPMagInt_Driver_ComputeAndSetErrcon( GPMagInt_Driver *This );

FQUALIFIER
void GPMagInt_Driver_ReSetParameters( GPMagInt_Driver *This,
                                      G4double new_safety);
FQUALIFIER
void GPMagInt_Driver_SetSafety( GPMagInt_Driver *This, G4double val);

FQUALIFIER
void GPMagInt_Driver_SetPgrow( GPMagInt_Driver *This, G4double  val);

FQUALIFIER
void GPMagInt_Driver_SetErrcon( GPMagInt_Driver *This, G4double val);

FQUALIFIER
void GPMagInt_Driver_RenewStepperAndAdjust( GPMagInt_Driver *This,
                                            GPClassicalRK4 *pItsStepper);
FQUALIFIER
void GPMagInt_Driver_SetChargeMomentumMass( GPMagInt_Driver *This,
                                            G4double particleCharge, // e+
                                            G4double MomentumXc,
                                            G4double Mass );
FQUALIFIER
const GPClassicalRK4* GPMagInt_Driver_GetStepper(GPMagInt_Driver *This);

FQUALIFIER
G4int GPMagInt_Driver_GetMaxNoSteps( GPMagInt_Driver *This );

FQUALIFIER
void GPMagInt_Driver_SetMaxNoSteps( GPMagInt_Driver *This, G4int val);

FQUALIFIER
void GPMagInt_Driver_GetDerivatives( GPMagInt_Driver *This,
                                     GPFieldTrack &y_curr,
                                     G4double     dydx[]);
FQUALIFIER
G4double GPMagInt_Driver_GetVerboseLevel( GPMagInt_Driver *This );

FQUALIFIER 
void GPMagInt_Driver_SetVerboseLevel( GPMagInt_Driver *This, G4int newLevel);

FQUALIFIER
G4double GPMagInt_Driver_GetSmallestFraction( GPMagInt_Driver *This );

}

#endif
