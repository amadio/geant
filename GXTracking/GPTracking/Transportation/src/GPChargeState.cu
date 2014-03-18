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
// $Id: G4FieldTrack.icc,v 1.21 2006-11-13 18:24:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

// Implementation methods for the embedded class G4ChargeState
//                                ----------------------------

#include "GPChargeState.h"
#include "GPConstants.h"

FQUALIFIER 
void GPChargeState_Constructor(GPChargeState *This,
			       G4double charge,
			       G4double magnetic_dipole_moment,  
			       G4double electric_dipole_moment,  
			       G4double magnetic_charge)
{
   This->fCharge= charge;
   This->fMagn_dipole= magnetic_dipole_moment;
   This->fElec_dipole= electric_dipole_moment;
   This->fMagneticCharge= magnetic_charge;    
}  

FQUALIFIER
void GPChargeState_Copy_Constructor(GPChargeState *This, GPChargeState& Right )
{ 
  This->fCharge= Right.fCharge;
  This->fMagn_dipole= Right.fMagn_dipole; 
  This->fElec_dipole= Right.fElec_dipole;
  This->fMagneticCharge= Right.fMagneticCharge;
}

FQUALIFIER 
void GPChargeState_SetChargeAndMoments(GPChargeState *This,
				       G4double charge, 
				       G4double magnetic_dipole_moment,  
				       G4double electric_dipole_moment,
				       G4double magnetic_charge )
   //  Revise the charge and potentially all moments.
   //   By default do not change mdm, edm, mag charge. 
{
   This->fCharge= charge;
   if( magnetic_dipole_moment < DBL_MAX) 
     This->fMagn_dipole= magnetic_dipole_moment;
   if( electric_dipole_moment < DBL_MAX) 
     This->fElec_dipole= electric_dipole_moment;
   if( magnetic_charge < DBL_MAX)        
     This->fMagneticCharge= magnetic_charge;    
}

FQUALIFIER
void GPChargeState_SetCharge(GPChargeState *This,
			     G4double charge)
{ 
  This->fCharge= charge; 
}

FQUALIFIER
G4double GPChargeState_GetCharge(GPChargeState *This) 
{ 
  return This->fCharge; 
}  

FQUALIFIER
G4double GPChargeState_GetMagneticDipoleMoment(GPChargeState *This) 
{ 
  return This->fMagn_dipole; 
} 

FQUALIFIER
G4double GPChargeState_ElectricDipoleMoment(GPChargeState *This) 
{ 
  return This->fElec_dipole; 
} 

FQUALIFIER
G4double GPChargeState_MagneticCharge(GPChargeState *This) 
{ 
  return This->fMagneticCharge; 
} 
