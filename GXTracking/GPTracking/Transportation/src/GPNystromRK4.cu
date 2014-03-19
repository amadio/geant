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
// $Id: GPNystromRK4.cc,v 1.9 2010-09-10 15:42:09 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// History:
// - Created:      I.Gavrilenko    15.05.2009   (as G4AtlasRK4)
// - Adaptations:  J. Apostolakis  May-Nov 2009
// -------------------------------------------------------------------

#include "GPNystromRK4.h"

//////////////////////////////////////////////////////////////////
// Constructor - with optional distance ( has default value)
//////////////////////////////////////////////////////////////////

FQUALIFIER
void GPNystromRK4_Constructor(GPNystromRK4 *This, 
			      GPEquationOfMotion* EqRhs, 
			      G4double distanceConstField)
{

  //  : G4MagIntegratorStepper(magEqRhs, 6),            // number of variables

  //EqRhs should be the one after GPEquationOfMotion_SetChargeMomentumMass is called

  GPNystromRK4_G4MagIntegratorStepper_Constructor(This,EqRhs,6,0);

  This->m_magdistance = distanceConstField;
  This->m_cof = 0.0;
  This->m_mom = 0.0;
  This->m_imom = 0.0;
  This->m_cachedMom = false; 

  This->m_fldPosition[0]  = This->m_iPoint[0] = This->m_fPoint[0] = This->m_mPoint[0] = 9.9999999e+99 ;
  This->m_fldPosition[1]  = This->m_iPoint[1] = This->m_fPoint[1] = This->m_mPoint[1] = 9.9999999e+99 ;
  This->m_fldPosition[2]  = This->m_iPoint[2] = This->m_fPoint[2] = This->m_mPoint[2] = 9.9999999e+99 ;
  This->m_fldPosition[3]  = -9.9999999e+99;
  This->m_lastField[0] = This->m_lastField[1] = This->m_lastField[2] = 0.0;

  This->m_magdistance2 = distanceConstField*distanceConstField;
}

////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////
/*
GPNystromRK4_Destructor()
{
}
*/
/////////////////////////////////////////////////////////////////////////////////
// Integration in one  step 
/////////////////////////////////////////////////////////////////////////////////

FQUALIFIER
void GPNystromRK4_Stepper(GPNystromRK4 *This,
			  const G4double P[],const G4double dPdS[],
			  G4double Step,G4double Po[],G4double Err[])
{
  G4double R[3] = {   P[0],   P[1] ,    P[2]};
  G4double A[3] = {dPdS[0], dPdS[1], dPdS[2]};

  This->m_iPoint[0]=R[0]; This->m_iPoint[1]=R[1]; This->m_iPoint[2]=R[2];

  const G4double one_sixth= 1./6.;
  G4double S  =     Step   ;
  G4double S5 =  .5*Step   ;
  G4double S4 = .25*Step   ;
  G4double S6 =     Step * one_sixth;   // Step / 6.;


  // John A  added, in order to emulate effect of call to changed/derived RHS
  // This->m_mom   = sqrt(P[3]*P[3]+P[4]*P[4]+P[5]*P[5]); 
  // This->m_imom  = 1./This->m_mom;
  // This->m_cof   = m_fEq->FCof()*This->m_imom;

  // Point 1
  //
  G4double K1[3] = { This->m_imom*dPdS[3], This->m_imom*dPdS[4], This->m_imom*dPdS[5] };
  
  // Point2
  //
  G4double p[4] = {R[0]+S5*(A[0]+S4*K1[0]),
		   R[1]+S5*(A[1]+S4*K1[1]),
		   R[2]+S5*(A[2]+S4*K1[2]),
		   P[7]                   }; 

  GPNystromRK4_getField(This, p);

  G4double A2[3] = {A[0]+S5*K1[0],A[1]+S5*K1[1],A[2]+S5*K1[2]};
  G4double K2[3] = {(A2[1]*This->m_lastField[2]-A2[2]*This->m_lastField[1])*This->m_cof,
		    (A2[2]*This->m_lastField[0]-A2[0]*This->m_lastField[2])*This->m_cof,
		    (A2[0]*This->m_lastField[1]-A2[1]*This->m_lastField[0])*This->m_cof};
 
  This->m_mPoint[0]=p[0]; This->m_mPoint[1]=p[1]; This->m_mPoint[2]=p[2];

  // Point 3 with the same magnetic field
  //
  G4double A3[3] = {A[0]+S5*K2[0],A[1]+S5*K2[1],A[2]+S5*K2[2]};
  G4double K3[3] = {(A3[1]*This->m_lastField[2]-A3[2]*This->m_lastField[1])*This->m_cof,
		    (A3[2]*This->m_lastField[0]-A3[0]*This->m_lastField[2])*This->m_cof,
		    (A3[0]*This->m_lastField[1]-A3[1]*This->m_lastField[0])*This->m_cof};
  
  // Point 4
  //
  p[0] = R[0]+S*(A[0]+S5*K3[0]);
  p[1] = R[1]+S*(A[1]+S5*K3[1]);
  p[2] = R[2]+S*(A[2]+S5*K3[2]);             

  GPNystromRK4_getField(This,p);
  
  G4double A4[3] = {A[0]+S*K3[0],A[1]+S*K3[1],A[2]+S*K3[2]};
  G4double K4[3] = {(A4[1]*This->m_lastField[2]-A4[2]*This->m_lastField[1])*This->m_cof,
		    (A4[2]*This->m_lastField[0]-A4[0]*This->m_lastField[2])*This->m_cof,
		    (A4[0]*This->m_lastField[1]-A4[1]*This->m_lastField[0])*This->m_cof};
  
  // New position
  //
  Po[0] = P[0]+S*(A[0]+S6*(K1[0]+K2[0]+K3[0]));
  Po[1] = P[1]+S*(A[1]+S6*(K1[1]+K2[1]+K3[1]));
  Po[2] = P[2]+S*(A[2]+S6*(K1[2]+K2[2]+K3[2]));

  This->m_fPoint[0]=Po[0]; This->m_fPoint[1]=Po[1]; This->m_fPoint[2]=Po[2];

  // New direction
  //
  Po[3] = A[0]+S6*(K1[0]+K4[0]+2.*(K2[0]+K3[0]));
  Po[4] = A[1]+S6*(K1[1]+K4[1]+2.*(K2[1]+K3[1]));
  Po[5] = A[2]+S6*(K1[2]+K4[2]+2.*(K2[2]+K3[2]));

  // Errors
  //
  Err[3] = S*fabs(K1[0]-K2[0]-K3[0]+K4[0]);
  Err[4] = S*fabs(K1[1]-K2[1]-K3[1]+K4[1]);
  Err[5] = S*fabs(K1[2]-K2[2]-K3[2]+K4[2]);
  Err[0] = S*Err[3]                       ;
  Err[1] = S*Err[4]                       ;
  Err[2] = S*Err[5]                       ;
  Err[3]*= This->m_mom                          ;
  Err[4]*= This->m_mom                          ;
  Err[5]*= This->m_mom                          ;

  // Normalize momentum
  //
  G4double normF = This->m_mom/sqrt(Po[3]*Po[3]+Po[4]*Po[4]+Po[5]*Po[5]);
  Po [3]*=normF; Po[4]*=normF; Po[5]*=normF; 

  // Pass Energy, time unchanged -- time is not integrated !!
  Po[6]=P[6]; Po[7]=P[7];
}


/////////////////////////////////////////////////////////////////////////////////
// Estimate the maximum distance from the curve to the chord
/////////////////////////////////////////////////////////////////////////////////

FQUALIFIER
G4double GPNystromRK4_DistChord(GPNystromRK4 *This) 
{
  G4double ax = This->m_fPoint[0]-This->m_iPoint[0];  
  G4double ay = This->m_fPoint[1]-This->m_iPoint[1];  
  G4double az = This->m_fPoint[2]-This->m_iPoint[2];
  G4double dx = This->m_mPoint[0]-This->m_iPoint[0]; 
  G4double dy = This->m_mPoint[1]-This->m_iPoint[1]; 
  G4double dz = This->m_mPoint[2]-This->m_iPoint[2];
  G4double d2 = (ax*ax+ay*ay+az*az)    ; 

  if(d2!=0.) {
    G4double  s = (ax*dx+ay*dy+az*dz)/d2;
    dx         -= (s*ax)                ;
    dy         -= (s*ay)                ;
    dz         -= (s*az)                ;
  }
  return sqrt(dx*dx+dy*dy+dz*dz);
}

/////////////////////////////////////////////////////////////////////////////////
// Derivatives calculation - caching the momentum value
/////////////////////////////////////////////////////////////////////////////////

FQUALIFIER
void GPNystromRK4_ComputeRightHandSide(GPNystromRK4 *This, 
				       const G4double P[], G4double dPdS[])
{
  G4double P4vec[4]= { P[0], P[1], P[2], P[7] }; // Time is P[7]
  GPNystromRK4_getField(This,P4vec);
  This->m_mom   = sqrt(P[3]*P[3]+P[4]*P[4]+P[5]*P[5])     ; 
  This->m_imom  = 1./This->m_mom                                ;
  This->m_cof   = GPEquationOfMotion_FCof(This->m_fEq)*This->m_imom                    ;
  This->m_cachedMom = true                                ; // Caching the value
  dPdS[0] = P[3]*This->m_imom                             ; // dx /ds
  dPdS[1] = P[4]*This->m_imom                             ; // dy /ds
  dPdS[2] = P[5]*This->m_imom                             ; // dz /ds
  dPdS[3] = This->m_cof*(P[4]*This->m_lastField[2]-P[5]*This->m_lastField[1]) ; // dPx/ds
  dPdS[4] = This->m_cof*(P[5]*This->m_lastField[0]-P[3]*This->m_lastField[2]) ; // dPy/ds
  dPdS[5] = This->m_cof*(P[3]*This->m_lastField[1]-P[4]*This->m_lastField[0]) ; // dPz/ds
}

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////
FQUALIFIER
void  GPNystromRK4_SetDistanceForConstantField(GPNystromRK4 *This, G4double length )
{
  This->m_magdistance=   length;
  This->m_magdistance2 = length*length;
}

FQUALIFIER
G4double  GPNystromRK4_GetDistanceForConstantField(GPNystromRK4 *This)
{
  return This->m_magdistance; 
}

/////////////////////////////////////////////////////////////////////////////////
// Get value of magnetic field while checking distance from last stored call
/////////////////////////////////////////////////////////////////////////////////

FQUALIFIER
void GPNystromRK4_getField (GPNystromRK4 *This, const G4double P[4])
{
  
  G4double dx = P[0]-This->m_fldPosition[0];
  G4double dy = P[1]-This->m_fldPosition[1];
  G4double dz = P[2]-This->m_fldPosition[2];

  if((dx*dx+dy*dy+dz*dz) > This->m_magdistance2)
  {
    This->m_fldPosition[0] = P[0];
    This->m_fldPosition[1] = P[1];
    This->m_fldPosition[2] = P[2];
    This->m_fldPosition[3] = P[3];   //  Generally it is P[7] - changed convention !!
    GPEquationOfMotion_GetFieldValue(This->m_fEq, This->m_fldPosition, This->m_lastField);
  }
}

//--------------------------------------------------------------
// class G4MagIntegratorStepper
//--------------------------------------------------------------  

// G4MagIntegratorStepper::G4MagIntegratorStepper
FQUALIFIER
void GPNystromRK4_G4MagIntegratorStepper_Constructor(GPNystromRK4 *This,
						     GPEquationOfMotion* Equation,
						     G4int       num_integration_vars,
						     G4int       num_state_vars)
{
  This->m_fEq = Equation;
  This->fNoIntegrationVariables = num_integration_vars;
  This->fNoStateVariables = num_state_vars;
}
