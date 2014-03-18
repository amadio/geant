#include "GXNystromRK4.h"

FQUALIFIER
void GXNystromRK4_Constructor(GXNystromRK4 *This, 
			      GXEquationOfMotion* EqRhs) 
{

  This->m_fEq = EqRhs;

  This->m_fldPosition[0]  = This->m_iPoint[0] = This->m_fPoint[0] = This->m_mPoint[0] = 9.9999999e+99 ;
  This->m_fldPosition[1]  = This->m_iPoint[1] = This->m_fPoint[1] = This->m_mPoint[1] = 9.9999999e+99 ;
  This->m_fldPosition[2]  = This->m_iPoint[2] = This->m_fPoint[2] = This->m_mPoint[2] = 9.9999999e+99 ;
  This->m_lastField[0] = This->m_lastField[1] = This->m_lastField[2] = 0.0;

}

/////////////////////////////////////////////////////////////////////////////////
// Integration in one  step 
/////////////////////////////////////////////////////////////////////////////////

FQUALIFIER
void GXNystromRK4_Stepper(GXNystromRK4 *This,
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
  G4double m_mom   = sqrt(P[3]*P[3]+P[4]*P[4]+P[5]*P[5]); 
  G4double m_imom  = 1./m_mom;
  G4double m_cof   = GXEquationOfMotion_GetCoeff(This->m_fEq)*m_imom;

  // Point 1
  //
  G4double K1[3] = { m_imom*dPdS[3], m_imom*dPdS[4], m_imom*dPdS[5] };
  
  // Point2
  //
  G4double p[3] = {R[0]+S5*(A[0]+S4*K1[0]),
		   R[1]+S5*(A[1]+S4*K1[1]),
		   R[2]+S5*(A[2]+S4*K1[2])}; 

  GXNystromRK4_GetField(This, p);

  G4double A2[3] = {A[0]+S5*K1[0],A[1]+S5*K1[1],A[2]+S5*K1[2]};
  G4double K2[3] = {(A2[1]*This->m_lastField[2]-A2[2]*This->m_lastField[1])*m_cof,
		    (A2[2]*This->m_lastField[0]-A2[0]*This->m_lastField[2])*m_cof,
		    (A2[0]*This->m_lastField[1]-A2[1]*This->m_lastField[0])*m_cof};
 
  This->m_mPoint[0]=p[0]; This->m_mPoint[1]=p[1]; This->m_mPoint[2]=p[2];

  // Point 3 with the same magnetic field
  //
  G4double A3[3] = {A[0]+S5*K2[0],A[1]+S5*K2[1],A[2]+S5*K2[2]};
  G4double K3[3] = {(A3[1]*This->m_lastField[2]-A3[2]*This->m_lastField[1])*m_cof,
		    (A3[2]*This->m_lastField[0]-A3[0]*This->m_lastField[2])*m_cof,
		    (A3[0]*This->m_lastField[1]-A3[1]*This->m_lastField[0])*m_cof};
  
  // Point 4
  //
  p[0] = R[0]+S*(A[0]+S5*K3[0]);
  p[1] = R[1]+S*(A[1]+S5*K3[1]);
  p[2] = R[2]+S*(A[2]+S5*K3[2]);             

  GXNystromRK4_GetField(This,p);
  
  G4double A4[3] = {A[0]+S*K3[0],A[1]+S*K3[1],A[2]+S*K3[2]};
  G4double K4[3] = {(A4[1]*This->m_lastField[2]-A4[2]*This->m_lastField[1])*m_cof,
		    (A4[2]*This->m_lastField[0]-A4[0]*This->m_lastField[2])*m_cof,
		    (A4[0]*This->m_lastField[1]-A4[1]*This->m_lastField[0])*m_cof};
  
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
  Err[3]*= m_mom                          ;
  Err[4]*= m_mom                          ;
  Err[5]*= m_mom                          ;

  // Normalize momentum
  //
  G4double normF = m_mom/sqrt(Po[3]*Po[3]+Po[4]*Po[4]+Po[5]*Po[5]);
  Po [3]*=normF; Po[4]*=normF; Po[5]*=normF; 

}


/////////////////////////////////////////////////////////////////////////////////
// Estimate the maximum distance from the curve to the chord
/////////////////////////////////////////////////////////////////////////////////

FQUALIFIER
G4double GXNystromRK4_DistChord(GXNystromRK4 *This) 
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
void GXNystromRK4_RightHandSide(GXNystromRK4 *This, 
				const G4double P[], G4double dPdS[])
{
  GXEquationOfMotion_RightHandSide(This->m_fEq, P, dPdS);
}

/////////////////////////////////////////////////////////////////////////////////
// Get value of magnetic field while checking distance from last stored call
/////////////////////////////////////////////////////////////////////////////////

FQUALIFIER
void GXNystromRK4_GetField (GXNystromRK4 *This, const G4double P[4])
{
  This->m_fldPosition[0] = P[0];
  This->m_fldPosition[1] = P[1];
  This->m_fldPosition[2] = P[2];
  GXEquationOfMotion_GetFieldValue(This->m_fEq, This->m_fldPosition, This->m_lastField);
}

FQUALIFIER
G4double GXNystromRK4_TruncatedError(GXNystromRK4 *This,
				     G4double hstep,
				     const G4double yarrout[],
				     const G4double yerr_vec[] )
{
  G4double dyerr = 0.0;
  G4double dyerr_pos_sq, dyerr_mom_rel_sq;  
  G4double dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

  vel_mag_sq = (yarrout[3]*yarrout[3] + yarrout[4]*yarrout[4] +
  		 yarrout[5]*yarrout[5] );
  inv_vel_mag_sq = 1.0 / vel_mag_sq; 

  dyerr_pos_sq = ( yerr_vec[0]*yerr_vec[0]+yerr_vec[1]*yerr_vec[1]+
          yerr_vec[2]*yerr_vec[2]);
  dyerr_mom_sq = ( yerr_vec[3]*yerr_vec[3]+yerr_vec[4]*yerr_vec[4]+
          yerr_vec[5]*yerr_vec[5]);
  dyerr_mom_rel_sq = dyerr_mom_sq * inv_vel_mag_sq;

  if( dyerr_pos_sq > ( dyerr_mom_rel_sq * hstep*hstep ) )
  {
    dyerr = sqrt(dyerr_pos_sq);
  }
  else
  {
    dyerr = sqrt(dyerr_mom_rel_sq) * hstep;
  }

  return dyerr;
}
