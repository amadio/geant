#include "GXCashKarpRKF45.h"
#include "GPUtils.h"

FQUALIFIER
void GXCashKarpRKF45_Constructor( GXCashKarpRKF45 *This, 
				  GXEquationOfMotion* EqRhs)
{
  This->fEquation_Rhs = EqRhs;
  This->fLastStepLength = 0.0;
}

FQUALIFIER
G4int GXCashKarpRKF45_IntegratorOrder() { return 4; }

//////////////////////////////////////////////////////////////////////
//
// Given values for n = 6 variables yIn[0,...,n-1] 
// known  at x, use the fifth-order Cash-Karp Runge-
// Kutta-Fehlberg-4-5 method to advance the solution over an interval
// Step and return the incremented variables as yOut[0,...,n-1]. Also
// return an estimate of the local truncation error yErr[] using the
// embedded 4th-order method. The user supplies routine
// RightHandSide(y,dydx), which returns derivatives dydx for y .

FQUALIFIER
void GXCashKarpRKF45_Stepper( GXCashKarpRKF45 *This, 
			      const G4double yInput[],
			      const G4double dydx[],
			      G4double Step,
			      G4double yOut[],
			      G4double yErr[])

{  
  G4int i;

  const G4double  b21 = 0.2 ,
    b31 = 3.0/40.0 , b32 = 9.0/40.0 ,
    b41 = 0.3 , b42 = -0.9 , b43 = 1.2 ,
    
    b51 = -11.0/54.0 , b52 = 2.5 , b53 = -70.0/27.0 ,
    b54 = 35.0/27.0 ,
    
    b61 = 1631.0/55296.0 , b62 =   175.0/512.0 ,
    b63 =  575.0/13824.0 , b64 = 44275.0/110592.0 ,
    b65 =  253.0/4096.0 ,
    
    c1 = 37.0/378.0 , c3 = 250.0/621.0 , c4 = 125.0/594.0 ,
    c6 = 512.0/1771.0 ,
    dc5 = -277.0/14336.0 ;
  
  const G4double dc1 = c1 - 2825.0/27648.0 ,  dc3 = c3 - 18575.0/48384.0 ,
    dc4 = c4 - 13525.0/55296.0 , dc6 = c6 - 0.25 ;
  
  // The number of variables to be integrated over
  const G4int numberOfVariables = 6;

  // scratch space moved from the head file
  G4double ak2[6], ak3[6], ak4[6], ak5[6], ak6[6];
  G4double yTemp[6],yIn[6];

  //  Saving yInput because yInput and yOut can be aliases for same array

  for(i=0;i<numberOfVariables;i++) {
    yIn[i]=yInput[i];
  }
  // RightHandSide(This,yIn, dydx) ;              // 1st Step

  for(i=0;i<numberOfVariables;i++) {
    yTemp[i] = yIn[i] + b21*Step*dydx[i] ;
  }
  GXCashKarpRKF45_RightHandSide(This, yTemp, ak2) ;              // 2nd Step

  for(i=0;i<numberOfVariables;i++) {
    yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
  }
  GXCashKarpRKF45_RightHandSide(This, yTemp, ak3) ;              // 3rd Step

  for(i=0;i<numberOfVariables;i++) {
    yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] 
					  + b43*ak3[i]) ;
  }
  GXCashKarpRKF45_RightHandSide(This, yTemp, ak4) ;              // 4th Step
  
  for(i=0;i<numberOfVariables;i++) {
    yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] 
		                  + b53*ak3[i] + b54*ak4[i]) ;
  }
  GXCashKarpRKF45_RightHandSide(This, yTemp, ak5) ;              // 5th Step

  for(i=0;i<numberOfVariables;i++) {
    yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] 
		   + b63*ak3[i] + b64*ak4[i] + b65*ak5[i]) ;
  }
  GXCashKarpRKF45_RightHandSide(This, yTemp, ak6) ;              // 6th Step

  for(i=0;i<numberOfVariables;i++) {
    // Accumulate increments with proper weights

    yOut[i] = yIn[i] + Step*(c1*dydx[i] + c3*ak3[i] 
			   + c4*ak4[i] + c6*ak6[i]) ;
    
    // Estimate error as difference between 4th and
    // 5th order methods
    
    yErr[i] = Step*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] +
		    dc5*ak5[i] + dc6*ak6[i]) ;
    
    // Store Input and Final values, for possible use in calculating chord
    This->fLastInitialVector[i] = yIn[i] ;
    This->fLastFinalVector[i]   = yOut[i];
    This->fLastDyDx[i]          = dydx[i];
  }
  // NormaliseTangentVector( yOut ); // Not wanted
  
  This->fLastStepLength =Step;

 return ;
}

FQUALIFIER 
void GXCashKarpRKF45_RightHandSide(GXCashKarpRKF45 *This,
				  const  G4double y[], 
				  G4double dydx[] )   
{
  GXEquationOfMotion_RightHandSide(This->fEquation_Rhs, y, dydx);
}

// G4MagErrorStepper::DistChord
FQUALIFIER
G4double GXCashKarpRKF45_DistChord(GXCashKarpRKF45 *This)
{
  // Estimate the maximum distance from the curve to the chord
  G4double distLine, distChord; 
  GPThreeVector initialPoint, finalPoint, midPoint;

  // Store last initial and final points (they will be overwritten 
  // in self-Stepper call!)
  initialPoint = GPThreeVector_create( This->fLastInitialVector[0], 
				       This->fLastInitialVector[1], 
				       This->fLastInitialVector[2]); 
  finalPoint   = GPThreeVector_create( This->fLastFinalVector[0],  
				       This->fLastFinalVector[1],  
				       This->fLastFinalVector[2]); 
  
  // Do half a step using StepNoErr

  GXCashKarpRKF45_Stepper(This, This->fLastInitialVector, 
			  This->fLastDyDx, 0.5 * This->fLastStepLength, 
			  This->fMidVector, This->fMidError );

  midPoint = GPThreeVector_create( This->fMidVector[0], 
				   This->fMidVector[1], 
				   This->fMidVector[2]);       

  // Use stored values of Initial and Endpoint + new Midpoint to evaluate
  //  distance of Chord

  if (GPThreeVector_nequal(initialPoint,finalPoint)) {
    distLine= GPLineSection_Distline( midPoint,initialPoint, 
				      finalPoint );
    distChord = distLine;
  }else{
    distChord = GPThreeVector_mag(GPThreeVector_sub(midPoint,
						    initialPoint));
  }
  return distChord;
}

FQUALIFIER
G4double GXCashKarpRKF45_TruncatedError(GXCashKarpRKF45 *This,
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
