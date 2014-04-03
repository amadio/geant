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
// $Id: G4FieldTrack.hh,v 1.21 2006-11-13 18:24:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4FieldTrack
//
// Class description:
//
// Data structure bringing together a magnetic track's state.
// (position, momentum direction & modulus, energy, spin, ... )
// Uses/abilities:
//  - does not maintain any relationship between its data (eg energy/momentum).
//  - for use in Runge-Kutta solver (in passing it the values right now).

// History
// - First version: Oct 14, 1996  John Apostolakis
// - Modified:      Oct 24, 1996  JA: Added dist_on_curve, deleted constructor
//                  Nov  5, 1998  JA: Added energy, momentum, TOF, spin &
//                                    several constructor, access, set methods
//                  May 10, 2006  JA: Added charge, "default" constructor
// -------------------------------------------------------------------

#ifndef GPFieldTrack_HH
#define GPFieldTrack_HH

#include "GPTypeDef.h"
#include "GPConstants.h"
#include "GPThreeVector.h"
#include "GPChargeState.h"

struct GPFieldTrack 
{
  G4double  SixVector[6];
  G4double  fDistanceAlongCurve;  // distance along curve of point
  G4double  fKineticEnergy;
  G4double  fRestMass_c2;
  G4double  fLabTimeOfFlight;
  G4double  fProperTimeOfFlight;
  GPThreeVector fSpin;
  GPThreeVector fMomentumDir;

  GPChargeState fChargeState;

};

const int GPFieldTrack_ncompSVEC = 12;

extern "C" {

FQUALIFIER
void GPFieldTrack_Constructor( GPFieldTrack *This, 
                               const GPThreeVector& pPosition, 
                               G4double       LaboratoryTimeOfFlight,
                               const GPThreeVector& pMomentumDirection,
                               G4double       kineticEnergy,
                               G4double       restMass_c2,
                               G4double       charge, 
                               const GPThreeVector& Spin,
                               G4double       magnetic_dipole_moment=0.0,
                               G4double       curve_length=0.0 );

FQUALIFIER
void GPFieldTrack_Constructor2( GPFieldTrack *This, 
                                const GPThreeVector& pPosition, 
                                const GPThreeVector& pMomentumDirection,    
                                G4double       curve_length, 
                                G4double       kineticEnergy,
                                const G4double       restMass_c2,
                                G4double       velocity=0.0,
                                G4double       pLaboratoryTimeOfFlight=0.0,
                                G4double       pProperTimeOfFlight=0.0,
                                const GPThreeVector* pSpin=0);

FQUALIFIER
void GPFieldTrack_SetChargeAndMoments(GPFieldTrack *This,
				      G4double charge, 
				      G4double magnetic_dipole_moment= DBL_MAX, 
				      G4double electric_dipole_moment= DBL_MAX, 
				      G4double magnetic_charge= DBL_MAX );

FQUALIFIER
void GPFieldTrack_InitialiseSpin(GPFieldTrack *This, 
				 const GPThreeVector& Spin );

FQUALIFIER
void GPFieldTrack_UpdateFourMomentum( GPFieldTrack *This,
				      G4double             kineticEnergy, 
				      const GPThreeVector& momentumDirection );

FQUALIFIER
void GPFieldTrack_UpdateState(GPFieldTrack *This, 
			      const GPThreeVector& position, 
			      G4double             laboratoryTimeOfFlight,
			      const GPThreeVector& momentumDirection,
			      G4double             kineticEnergy);
                              
FQUALIFIER
GPFieldTrack& GPFieldTrack_SetCurvePnt(GPFieldTrack *This, 
				       const GPThreeVector& pPosition, 
				       const GPThreeVector& pMomentum,  
				       G4double       s_curve );

FQUALIFIER
GPThreeVector  GPFieldTrack_GetMomentum(const GPFieldTrack *This);

FQUALIFIER
GPThreeVector  GPFieldTrack_GetPosition(const GPFieldTrack *This);

FQUALIFIER
GPThreeVector GPFieldTrack_GetMomentumDir(const GPFieldTrack *This);

FQUALIFIER
GPThreeVector  GPFieldTrack_GetMomentumDirection(const GPFieldTrack *This);

FQUALIFIER
G4double       GPFieldTrack_GetCurveLength(const GPFieldTrack *This);

FQUALIFIER
GPThreeVector  GPFieldTrack_GetSpin(GPFieldTrack *This);

FQUALIFIER
G4double       GPFieldTrack_GetLabTimeOfFlight(GPFieldTrack *This);

FQUALIFIER
G4double GPFieldTrack_GetKineticEnergy(GPFieldTrack *This);

FQUALIFIER
G4double GPFieldTrack_GetCharge(GPFieldTrack *This);

FQUALIFIER
void GPFieldTrack_SetPosition(GPFieldTrack *This, GPThreeVector pPosition) ;

FQUALIFIER
void GPFieldTrack_SetMomentum(GPFieldTrack *This, GPThreeVector pMomentum );

FQUALIFIER
void GPFieldTrack_SetMomentumDir(GPFieldTrack *This, GPThreeVector newMomDir);

FQUALIFIER
void GPFieldTrack_SetCurveLength(GPFieldTrack *This, G4double nCurve_s);

FQUALIFIER
void GPFieldTrack_SetKineticEnergy(GPFieldTrack *This, G4double newKinEnergy);

FQUALIFIER
void GPFieldTrack_SetSpin(GPFieldTrack *This, GPThreeVector nSpin);

FQUALIFIER
void GPFieldTrack_SetLabTimeOfFlight(GPFieldTrack *This, G4double nTOF) ;

FQUALIFIER
void GPFieldTrack_SetProperTimeOfFlight(GPFieldTrack *This, G4double nTOF);

FQUALIFIER
void GPFieldTrack_DumpToArray(GPFieldTrack *This,
			      G4double valArr[GPFieldTrack_ncompSVEC] );

FQUALIFIER
void GPFieldTrack_LoadFromArray(GPFieldTrack *This,
			       const G4double valArrIn[GPFieldTrack_ncompSVEC], 
				G4int noVarsIntegrated);

FQUALIFIER
GPChargeState* GetChargeState(GPFieldTrack *This) ;

}

#endif 
