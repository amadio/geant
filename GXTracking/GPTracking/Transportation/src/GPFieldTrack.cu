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
// $Id: G4FieldTrack.cc,v 1.15 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

#include "GPFieldTrack.h"

FQUALIFIER
void GPFieldTrack_Constructor( GPFieldTrack *This, 
			       const GPThreeVector& pPosition, 
			       G4double       LaboratoryTimeOfFlight,
			       const GPThreeVector& pMomentumDirection,
			       G4double       kineticEnergy,
			       G4double       restMass_c2,
			       G4double       charge, 
			       const GPThreeVector& Spin,
			       G4double       magnetic_dipole_moment,
			       G4double       curve_length )
{
  This->fKineticEnergy = kineticEnergy;
  This->fRestMass_c2 = restMass_c2;
  This->fLabTimeOfFlight = LaboratoryTimeOfFlight; 
  This->fProperTimeOfFlight = 0.;

  //  fChargeState(  charge, magnetic_dipole_moment ) 
  GPChargeState_Constructor(&This->fChargeState, charge, magnetic_dipole_moment,0,0); 

  G4double momentum  = sqrt(kineticEnergy*(kineticEnergy+2.0*restMass_c2));
  GPThreeVector pMomentum= GPThreeVector_mult(pMomentumDirection,momentum); 
  GPFieldTrack_SetCurvePnt(This, pPosition, pMomentum, curve_length );
  This->fMomentumDir=pMomentumDirection; 
  GPFieldTrack_InitialiseSpin(This, Spin ); 
}

FQUALIFIER
void GPFieldTrack_Constructor2( GPFieldTrack *This, 
				const GPThreeVector& pPosition, 
				const GPThreeVector& pMomentumDirection,    
				G4double       curve_length, 
				G4double       kineticEnergy,
				const G4double       restMass_c2,
				G4double       velocity,
				G4double       pLaboratoryTimeOfFlight,
				G4double       pProperTimeOfFlight,
				const GPThreeVector* pSpin)
{
  This->fKineticEnergy = kineticEnergy;
  This->fRestMass_c2 = restMass_c2;
  This->fLabTimeOfFlight = pLaboratoryTimeOfFlight; 
  This->fProperTimeOfFlight = pProperTimeOfFlight;
  //  This->fChargeState = DBL_MAX;
  GPChargeState_Constructor(&This->fChargeState,DBL_MAX,0,0,0); 

  G4double momentum  = sqrt(kineticEnergy*(kineticEnergy+2.0*restMass_c2));
  GPThreeVector pMomentum= GPThreeVector_mult(pMomentumDirection,momentum); 

  GPFieldTrack_SetCurvePnt(This, pPosition, pMomentum, curve_length );

  This->fMomentumDir=pMomentumDirection; 

  GPThreeVector Spin; 
  if( !pSpin ) Spin = GPThreeVector_create(0.,0.,0.); 
  else         Spin = *pSpin;

  GPFieldTrack_InitialiseSpin(This, Spin ); 
}

//FQUALIFIER
//void GPFieldTrack_Copy_Constructor(GPFieldTrack *this,
//				   GPFieldTrack&  rStVec)
//{
  /*
  This->fDistanceAlongCurve = rStVec.fDistanceAlongCurve;
  This->fKineticEnergy = rStVec.fKineticEnergy;
  This->fRestMass_c2 = rStVec.fRestMass_c2;
  This->fLabTimeOfFlight = rStVec.fLabTimeOfFlight; 
  This->fProperTimeOfFlight = rStVec.fProperTimeOfFlight; 
   // fMomentumModulus( rStVec.fMomentumModulus ),
  This->fSpin = rStVec.fSpin; 
  This->fMomentumDir = rStVec.fMomentumDir;

  This->fChargeState = rStVec.fChargeState;

  This->SixVector[0]= rStVec.SixVector[0];
  This->SixVector[1]= rStVec.SixVector[1];
  This->SixVector[2]= rStVec.SixVector[2];
  This->SixVector[3]= rStVec.SixVector[3];
  This->SixVector[4]= rStVec.SixVector[4];
  This->SixVector[5]= rStVec.SixVector[5];
  */

  // fpChargeState= new G4ChargeState( *rStVec.fpChargeState );
  // Can share charge state only when using handles etc
  //   fpChargeState = rStVec.fpChargeState;  
//}

FQUALIFIER
void GPFieldTrack_SetChargeAndMoments(GPFieldTrack *This,
				      G4double charge, 
				      G4double magnetic_dipole_moment, 
				      G4double electric_dipole_moment, 
				      G4double magnetic_charge)
{
  GPChargeState_SetChargeAndMoments(&This->fChargeState, 
				    charge,  magnetic_dipole_moment, 
				    electric_dipole_moment,  magnetic_charge ); 
  
}

FQUALIFIER
void GPFieldTrack_InitialiseSpin(GPFieldTrack *This, 
				 const GPThreeVector& Spin )
{
  This->fSpin = Spin;
}

FQUALIFIER
void GPFieldTrack_UpdateFourMomentum( GPFieldTrack *This,
				  G4double             kineticEnergy, 
                                  const GPThreeVector& momentumDirection )
{
  G4double momentum_mag  = sqrt(kineticEnergy*kineticEnergy
                            +2.0*This->fRestMass_c2*kineticEnergy);
  GPThreeVector momentumVector= GPThreeVector_mult(momentumDirection,
						   momentum_mag); 

  GPFieldTrack_SetMomentum(This, momentumVector ); 
  This->fMomentumDir=   momentumDirection;
  This->fKineticEnergy= kineticEnergy;
}

FQUALIFIER
void GPFieldTrack_UpdateState(GPFieldTrack *This, 
			      const GPThreeVector& position, 
			      G4double             laboratoryTimeOfFlight,
			      const GPThreeVector& momentumDirection,
			      G4double             kineticEnergy )
{ 
  // SetCurvePnt( position, momentumVector, s_curve=0.0);     
  GPFieldTrack_SetPosition(This, position); 
  This->fLabTimeOfFlight= laboratoryTimeOfFlight;
  This->fDistanceAlongCurve= 0.0;

  GPFieldTrack_UpdateFourMomentum(This, kineticEnergy, momentumDirection); 
}

FQUALIFIER
GPFieldTrack& GPFieldTrack_SetCurvePnt(GPFieldTrack *This, 
				       const GPThreeVector& pPosition, 
				       const GPThreeVector& pMomentum,  
				       G4double       s_curve )
{
  This->SixVector[0] = GPThreeVector_x(pPosition); 
  This->SixVector[1] = GPThreeVector_y(pPosition); 
  This->SixVector[2] = GPThreeVector_z(pPosition); 

  This->SixVector[3] = GPThreeVector_x(pMomentum); 
  This->SixVector[4] = GPThreeVector_y(pMomentum); 
  This->SixVector[5] = GPThreeVector_z(pMomentum); 

  This->fMomentumDir = GPThreeVector_unit(pMomentum);

  This->fDistanceAlongCurve= s_curve;

  return *This;
} 

// Assignment operator

FQUALIFIER
GPThreeVector  GPFieldTrack_GetMomentum(const GPFieldTrack *This)
{
  return GPThreeVector_create( This->SixVector[3], 
			       This->SixVector[4], 
			       This->SixVector[5]);

}   

FQUALIFIER
GPThreeVector  GPFieldTrack_GetPosition(const GPFieldTrack *This)
{
  return GPThreeVector_create( This->SixVector[0], 
			       This->SixVector[1], 
			       This->SixVector[2]);
}

FQUALIFIER
GPThreeVector GPFieldTrack_GetMomentumDir(const GPFieldTrack *This)
{
  return This->fMomentumDir;
}

FQUALIFIER
GPThreeVector  GPFieldTrack_GetMomentumDirection(const GPFieldTrack *This)
{
  return This->fMomentumDir;
}

FQUALIFIER
G4double       GPFieldTrack_GetCurveLength(const GPFieldTrack *This)
{
  return This->fDistanceAlongCurve;  

}

FQUALIFIER
GPThreeVector  GPFieldTrack_GetSpin(GPFieldTrack *This)
{
  return This->fSpin;

}

FQUALIFIER
G4double       GPFieldTrack_GetLabTimeOfFlight(GPFieldTrack *This)
{
  return This->fLabTimeOfFlight;
}

FQUALIFIER
G4double       GPFieldTrack_GetProperTimeOfFlight(GPFieldTrack *This)
{
  return This->fProperTimeOfFlight;
}

FQUALIFIER
G4double GPFieldTrack_GetKineticEnergy(GPFieldTrack *This)
{
  return This->fKineticEnergy;
}

FQUALIFIER
G4double GPFieldTrack_GetCharge(GPFieldTrack *This)
{
  return GPChargeState_GetCharge(&This->fChargeState);
}
       
// Accessors.

FQUALIFIER
void GPFieldTrack_SetPosition(GPFieldTrack *This, GPThreeVector pPosition) 
{
   This->SixVector[0] = GPThreeVector_x(pPosition); 
   This->SixVector[1] = GPThreeVector_y(pPosition); 
   This->SixVector[2] = GPThreeVector_z(pPosition); 
} 

FQUALIFIER
void GPFieldTrack_SetMomentum(GPFieldTrack *This, GPThreeVector pMomentum )
{
  This->SixVector[3] = GPThreeVector_x(pMomentum); 
  This->SixVector[4] = GPThreeVector_y(pMomentum); 
  This->SixVector[5] = GPThreeVector_z(pMomentum); 
  This->fMomentumDir = GPThreeVector_unit(pMomentum); 
}

FQUALIFIER
void GPFieldTrack_SetMomentumDir(GPFieldTrack *This, GPThreeVector newMomDir)
{
   This->fMomentumDir = newMomDir;
}

FQUALIFIER
void GPFieldTrack_SetCurveLength(GPFieldTrack *This, G4double nCurve_s)
{
  This->fDistanceAlongCurve = nCurve_s;  
}

FQUALIFIER
void GPFieldTrack_SetKineticEnergy(GPFieldTrack *This, G4double newKinEnergy)
{
  This->fKineticEnergy=newKinEnergy;
}

FQUALIFIER
void GPFieldTrack_SetSpin(GPFieldTrack *This, GPThreeVector nSpin)
{
  This->fSpin=nSpin;
}

FQUALIFIER
void GPFieldTrack_SetLabTimeOfFlight(GPFieldTrack *This, G4double nTOF) 
{
  This->fLabTimeOfFlight=nTOF;
}

FQUALIFIER
void GPFieldTrack_SetProperTimeOfFlight(GPFieldTrack *This, G4double nTOF)
{
  This->fProperTimeOfFlight=nTOF;
}

FQUALIFIER
void GPFieldTrack_DumpToArray(GPFieldTrack *This,
			      G4double valArr[GPFieldTrack_ncompSVEC] )
{
  valArr[0]=This->SixVector[0];
  valArr[1]=This->SixVector[1];
  valArr[2]=This->SixVector[2];
  valArr[3]=This->SixVector[3];
  valArr[4]=This->SixVector[4];
  valArr[5]=This->SixVector[5];

  GPThreeVector Momentum = 
    GPThreeVector_create(valArr[3],valArr[4],valArr[5]);

  // The following components may or may not be integrated.
  valArr[6]= This->fKineticEnergy; 
  valArr[7]=This->fLabTimeOfFlight;
  valArr[8]=This->fProperTimeOfFlight;
  valArr[9]= GPThreeVector_x(This->fSpin);
  valArr[10]=GPThreeVector_y(This->fSpin);
  valArr[11]=GPThreeVector_z(This->fSpin);
}

FQUALIFIER
void GPFieldTrack_LoadFromArray(GPFieldTrack *This,
			       const G4double valArrIn[GPFieldTrack_ncompSVEC], 
			       G4int noVarsIntegrated)
{
  G4int i;

  // Fill the variables not integrated with zero -- so it's clear !!
  G4double valArr[GPFieldTrack_ncompSVEC];
  for( i=0; i<noVarsIntegrated; i++){
     valArr[i]= valArrIn[i];
  }
  for( i=noVarsIntegrated; i<GPFieldTrack_ncompSVEC; i++) {
     valArr[i]= 0.0; 
  }

  This->SixVector[0]=valArr[0];
  This->SixVector[1]=valArr[1];
  This->SixVector[2]=valArr[2];
  This->SixVector[3]=valArr[3];
  This->SixVector[4]=valArr[4];
  This->SixVector[5]=valArr[5];

  GPThreeVector Momentum = 
    GPThreeVector_create(valArr[3],valArr[4],valArr[5]);

  G4double momentum_square= GPThreeVector_mag2(Momentum);
  This->fMomentumDir= GPThreeVector_unit(Momentum);

  This->fKineticEnergy = momentum_square / 
    (sqrt(momentum_square+This->fRestMass_c2*This->fRestMass_c2)
     + This->fRestMass_c2 ); 

  This->fLabTimeOfFlight=valArr[7];
  This->fProperTimeOfFlight=valArr[8];
  This->fSpin=GPThreeVector_create(valArr[9],valArr[10],valArr[11]);

}

FQUALIFIER
GPChargeState* GetChargeState(GPFieldTrack *This) 
{ 
  return &This->fChargeState; 
}
