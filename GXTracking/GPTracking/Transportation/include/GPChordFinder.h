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
// $Id: G4ChordFinder.hh,v 1.21 2008-10-29 14:17:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4ChordFinder
//
// class description:
//
// A class that provides RK integration of motion ODE  (as does g4magtr)
// and also has a method that returns an Approximate point on the curve 
// near to a (chord) point.

// History:
// - 25.02.97 - John Apostolakis - Design and implementation 
// -------------------------------------------------------------------

#ifndef GPCHORDFINDER_HH
#define GPCHORDFINDER_HH

#include "GPMagIntegratorDriver.h"
#include "GPClassicalRK4.h"
#include "GPEquationOfMotion.h"
#include "GPFieldTrack.h"
#include "GPChordFinder.h"

struct GPChordFinder
{ 
  // Constants
  G4double fDefaultDeltaChord;  // SET in G4ChordFinder.cc = 0.25 mm
  
  //  PARAMETERS 
  //  ---------------------
  G4double  fDeltaChord;               //  Maximum miss distance 
  //    Internal parameters
  G4double  fFirstFraction, fFractionLast, fFractionNextEstimate;
  G4double  fMultipleRadius; 
  G4int     fStatsVerbose;  // if > 0, print Statistics in destructor
  
  //  DEPENDENT Objects
  //  ---------------------
  GPMagInt_Driver*        fIntgrDriver;
  //  G4MagIntegratorStepper* fDriversStepper; 
  GPClassicalRK4* fDriversStepper; 
  G4bool                  fAllocatedStepper; // Bookkeeping of dependent object
  GPEquationOfMotion*     fEquation; 
  
  //  STATE information
  //  --------------------
  G4double    fLastStepEstimate_Unconstrained;
  //  State information for efficiency
  
  // For Statistics
  // -- G4int   fNoTrials, fNoCalls;
  G4int fTotalNoTrials_FNC, fNoCalls_FNC, fmaxTrials_FNC;//fnoTimesMaxTrFNC; 
};

extern "C" {

FQUALIFIER
void GPChordFinder_Constructor(GPChordFinder *This,
			       GPMagInt_Driver* pIntegrationDriver);

FQUALIFIER
void GPChordFinder_Constructor2( GPChordFinder *This,
				 GPMagneticField*        theMagField,
				 G4double                stepMinimum, 
				 GPClassicalRK4* pItsStepper );

FQUALIFIER
void   
GPChordFinder_SetFractions_Last_Next( GPChordFinder *This, 
				      G4double fractLast, G4double fractNext );

FQUALIFIER
G4double 
GPChordFinder_AdvanceChordLimited( GPChordFinder *This,
				   GPFieldTrack& yCurrent,
				   G4double      stepMax,
				   G4double      epsStep,
				   const GPThreeVector latestSafetyOrigin,
				   G4double       latestSafetyRadius );

FQUALIFIER
G4double
GPChordFinder_FindNextChord( GPChordFinder *This, 
			     const  GPFieldTrack& yStart,
			     G4double     stepMax,
			     GPFieldTrack&   yEnd, // Endpoint
			     G4double&   dyErrPos, // Error of endpoint
			     G4double    epsStep,
			     G4double*  pStepForAccuracy, 
			     const  GPThreeVector  latestSafetyOrigin,
			     G4double       latestSafetyRadius );

FQUALIFIER
G4double GPChordFinder_NewStep( GPChordFinder *This,
				G4double  stepTrialOld, 
                                G4double  dChordStep, // Curr. dchord achieved
                                G4double& stepEstimate_Unconstrained ) ;

FQUALIFIER
GPFieldTrack
GPChordFinder_ApproxCurvePointS( GPChordFinder *This, 
				 const GPFieldTrack&  CurveA_PointVelocity, 
				 const GPFieldTrack&  CurveB_PointVelocity, 
				 const GPFieldTrack&  ApproxCurveV,
				 const GPThreeVector& CurrentE_Point,
				 const GPThreeVector& CurrentF_Point,
				 const GPThreeVector& PointG,
				 G4bool first, G4double eps_step);

FQUALIFIER
GPFieldTrack 
GPChordFinder_ApproxCurvePointV( GPChordFinder *This, 
				 const GPFieldTrack& CurveA_PointVelocity, 
				 const GPFieldTrack& CurveB_PointVelocity, 
				 const GPThreeVector& CurrentE_Point,
				 G4double eps_step);

FQUALIFIER 
void GPChordFinder_SetIntegrationDriver(GPChordFinder *This,
					GPMagInt_Driver* IntegrationDriver);

FQUALIFIER
GPMagInt_Driver* GPChordFinder_GetIntegrationDriver(GPChordFinder *This);

FQUALIFIER
G4bool GPChordFinder_AcceptableMissDist(GPChordFinder *This,
					G4double dChordStep);

FQUALIFIER
void GPChordFinder_SetChargeMomentumMass(GPChordFinder *This,
					 G4double pCharge, // in e+ units
					 G4double pMomentum,
					 G4double pMass);

FQUALIFIER
G4double GPChordFinder_GetDeltaChord(GPChordFinder *This);

FQUALIFIER
void GPChordFinder_SetDeltaChord(GPChordFinder *This,G4double newval);

FQUALIFIER 
G4double GPChordFinder_GetLastStepEstimateUnc(GPChordFinder *This);

FQUALIFIER 
void GPChordFinder_SetLastStepEstimateUnc(GPChordFinder *This,
					  G4double stepEst );

FQUALIFIER
void GPChordFinder_ResetStepEstimate(GPChordFinder *This);

FQUALIFIER 
G4int GPChordFinder_GetNoCalls(GPChordFinder *This);

FQUALIFIER 
G4int GPChordFinder_GetNoTrials(GPChordFinder *This);

FQUALIFIER 
G4int GPChordFinder_GetNoMaxTrials(GPChordFinder *This);

FQUALIFIER 
G4double GPChordFinder_GetFirstFraction(GPChordFinder *This);

FQUALIFIER 
G4double GPChordFinder_GetFractionLast(GPChordFinder *This);

FQUALIFIER 
G4double GPChordFinder_GetFractionNextEstimate(GPChordFinder *This);

FQUALIFIER 
G4double GPChordFinder_GetMultipleRadius(GPChordFinder *This);

FQUALIFIER 
void GPChordFinder_SetFirstFraction(GPChordFinder *This,
				    G4double val);

FQUALIFIER 
G4int GPChordFinder_SetVerbose(GPChordFinder *This, G4int newvalue );

FQUALIFIER  
void GPChordFinder_AccumulateStatistics(GPChordFinder *This, 
					G4int noTrials );

FQUALIFIER 
G4double GPChordFinder_InvParabolic ( const G4double xa, 
				      const G4double ya,
				      const G4double xb, 
				      const G4double yb,
				      const G4double xc, 
				      const G4double yc );

}

#endif
