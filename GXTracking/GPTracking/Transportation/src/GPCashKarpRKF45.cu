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
// $Id: G4CashKarpRKF45.cc,v 1.16 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// The Cash-Karp Runge-Kutta-Fehlberg 4/5 method is an embedded fourth
// order method (giving fifth-order accuracy) for the solution of an ODE.
// Two different fourth order estimates are calculated; their difference
// gives an error estimate. [ref. Numerical Recipes in C, 2nd Edition]
// It is used to integrate the equations of the motion of a particle 
// in a magnetic field.
//
//  [ref. Numerical Recipes in C, 2nd Edition]
//
// -------------------------------------------------------------------

#include "GPCashKarpRKF45.h"
#include "GPUtils.h"

//--------------------------------------------------------------
// class hierarchy
//
// class G4ClassicalRK4 : public G4MagErrorStepper 
// class G4MagErrorStepper : public G4MagIntegratorStepper
// class G4MagIntegratorStepper: (stepper abstract base class). 
//--------------------------------------------------------------

//-------------------------------------------------
// class G4ClassicalRK4 : public G4MagErrorStepper 
//-------------------------------------------------

// G4ClassicalRK4::G4ClassicalRK4
FQUALIFIER
void GPCashKarpRKF45_Constructor( GPCashKarpRKF45 *This, 
				  GPEquationOfMotion* EqRhs,
				  G4int numberOfVariables,
				  G4bool primary)
{
  GPCashKarpRKF45_G4MagIntegratorStepper_Constructor(This,
		  EqRhs,numberOfVariables,0);

  This->fLastStepLength = 0.0;
}

FQUALIFIER
G4int GPCashKarpRKF45_IntegratorOrder() { return 4; }


// G4CashKarpRKF45::StepWithEst
FQUALIFIER
void GPCashKarpRKF45_StepWithEst( const G4double*,
				 const G4double*,
				 G4double,
				 G4double*,
				 G4double&,
				 G4double&,
				 const G4double*,
				 G4double*  ) 
{
  ;
  //  G4Exception("G4CashKarpRKF45::StepWithEst()", "GeomField0001",
  //              FatalException, "Method no longer used.");
  //  return ;

} 

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
void GPCashKarpRKF45_Stepper( GPCashKarpRKF45 *This, 
			      const G4double yInput[],
			      const G4double dydx[],
			      G4double Step,
			      G4double yOut[],
			      G4double yErr[])

{  
  // const G4int nvar = 6 ;
  // const G4double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
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
  
  // Initialise time to t0, needed when it is not updated by the integration.
  //        [ Note: Only for time dependent fields (usually electric) 
  //                  is it neccessary to integrate the time.] 
  yOut[7] = This->yTemp[7]   = This->yIn[7]; 
  
  const G4int numberOfVariables = GPCashKarpRKF45_GetNumberOfVariables(This) ;

  // The number of variables to be integrated over

  //  Saving yInput because yInput and yOut can be aliases for same array

  for(i=0;i<numberOfVariables;i++) {
    This->yIn[i]=yInput[i];
  }
  // RightHandSide(This,This->yIn, dydx) ;              // 1st Step

  for(i=0;i<numberOfVariables;i++) {
    This->yTemp[i] = This->yIn[i] + b21*Step*dydx[i] ;
  }
  GPCashKarpRKF45_RightHandSide(This, This->yTemp, This->ak2) ;              // 2nd Step

  for(i=0;i<numberOfVariables;i++) {
    This->yTemp[i] = This->yIn[i] + Step*(b31*dydx[i] + b32*This->ak2[i]) ;
  }
  GPCashKarpRKF45_RightHandSide(This, This->yTemp, This->ak3) ;              // 3rd Step

  for(i=0;i<numberOfVariables;i++) {
    This->yTemp[i] = This->yIn[i] + Step*(b41*dydx[i] + b42*This->ak2[i] 
					  + b43*This->ak3[i]) ;
  }
  GPCashKarpRKF45_RightHandSide(This, This->yTemp, This->ak4) ;              // 4th Step
  
  for(i=0;i<numberOfVariables;i++) {
    This->yTemp[i] = This->yIn[i] + Step*(b51*dydx[i] + b52*This->ak2[i] 
		                  + b53*This->ak3[i] + b54*This->ak4[i]) ;
  }
  GPCashKarpRKF45_RightHandSide(This, This->yTemp, This->ak5) ;              // 5th Step

  for(i=0;i<numberOfVariables;i++) {
    This->yTemp[i] = This->yIn[i] + Step*(b61*dydx[i] + b62*This->ak2[i] 
		   + b63*This->ak3[i] + b64*This->ak4[i] + b65*This->ak5[i]) ;
  }
  GPCashKarpRKF45_RightHandSide(This, This->yTemp, This->ak6) ;              // 6th Step

  for(i=0;i<numberOfVariables;i++) {
    // Accumulate increments with proper weights

    yOut[i] = This->yIn[i] + Step*(c1*dydx[i] + c3*This->ak3[i] 
			   + c4*This->ak4[i] + c6*This->ak6[i]) ;
    
    // Estimate error as difference between 4th and
    // 5th order methods
    
    yErr[i] = Step*(dc1*dydx[i] + dc3*This->ak3[i] + dc4*This->ak4[i] +
		    dc5*This->ak5[i] + dc6*This->ak6[i]) ;
    
    // Store Input and Final values, for possible use in calculating chord
    This->fLastInitialVector[i] = This->yIn[i] ;
    This->fLastFinalVector[i]   = yOut[i];
    This->fLastDyDx[i]          = dydx[i];
  }
  // NormaliseTangentVector( yOut ); // Not wanted
  
  This->fLastStepLength =Step;

 return ;
}

// G4MagErrorStepper::DistChord
FQUALIFIER
G4double GPCashKarpRKF45_DistChord(GPCashKarpRKF45 *This)
{
  // Estimate the maximum distance from the curve to the chord
  //
  //  We estimate this using the distance of the midpoint to 
  //  chord (the line between 
  // 
  //  Method below is good only for angle deviations < 2 pi, 
  //   This restriction should not a problem for the Runge cutta methods, 
  //   which generally cannot integrate accurately for large angle deviations.
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

  GPCashKarpRKF45_Stepper(This, This->fLastInitialVector, 
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

//--------------------------------------------------------------
// class G4MagIntegratorStepper
//--------------------------------------------------------------  

// G4MagIntegratorStepper::G4MagIntegratorStepper
FQUALIFIER
void GPCashKarpRKF45_G4MagIntegratorStepper_Constructor(GPCashKarpRKF45 *This,
				GPEquationOfMotion* Equation,
				G4int       num_integration_vars,
				G4int       num_state_vars)
{
  This->fEquation_Rhs = Equation;
  This->fNoIntegrationVariables = num_integration_vars;
  This->fNoStateVariables = num_state_vars;
}

// G4MagIntegratorStepper::ComputeRightHandSide
FQUALIFIER
void GPCashKarpRKF45_ComputeRightHandSide( GPCashKarpRKF45 *This,
					  const G4double y[], 
					  G4double dydx[] ) 
{
  GPCashKarpRKF45_RightHandSide(This, y, dydx );
}

// G4MagIntegratorStepper::GetEquationOfMotion
FQUALIFIER
GPEquationOfMotion* GPCashKarpRKF45_GetEquationOfMotion(GPCashKarpRKF45 *This)
{
  return  This->fEquation_Rhs;
} 

// G4MagIntegratorStepper::SetEquationOfMotion
FQUALIFIER
void GPCashKarpRKF45_SetEquationOfMotion(GPCashKarpRKF45 *This,
					GPEquationOfMotion* newEquation)
{
  if( newEquation != 0 )
  {
    This->fEquation_Rhs= newEquation;
  }
} 

// G4MagIntegratorStepper::GetNumberOfVariables
FQUALIFIER
G4int GPCashKarpRKF45_GetNumberOfVariables(GPCashKarpRKF45 *This) 
{
  return This->fNoIntegrationVariables;
}

// G4MagIntegratorStepper::GetNumberOfStateVariables
FQUALIFIER
G4int GPCashKarpRKF45_GetNumberOfStateVariables(GPCashKarpRKF45 *This) 
{
  return This->fNoStateVariables;
}

// G4MagIntegratorStepper::RightHandSide
FQUALIFIER 
void GPCashKarpRKF45_RightHandSide(GPCashKarpRKF45 *This,
				  const  G4double y[], 
				  G4double dydx[] )   
{
  GPEquationOfMotion_RightHandSide(This->fEquation_Rhs, y, dydx);
}

// G4MagIntegratorStepper::NormaliseTangentVector
FQUALIFIER 
void GPCashKarpRKF45_NormaliseTangentVector(GPCashKarpRKF45 *This, 
					   G4double vec[6] )
{
  G4double drds2 = vec[3]*vec[3]+vec[4]*vec[4]+vec[5]*vec[5];

  if( fabs(drds2 - 1.0) > 1.e-14 )
  {
    G4double normx = 1.0 / sqrt(drds2);
    for(G4int i=3;i<6;i++) { vec[i] *= normx; }
  }
}

// G4MagIntegratorStepper::NormalisePolarizationVector
FQUALIFIER 
void GPCashKarpRKF45_NormalisePolarizationVector(GPCashKarpRKF45 *This, 
						G4double vec[12] )
{
  G4double drds2 = vec[9]*vec[9]+vec[10]*vec[10]+vec[11]*vec[11];
  
  if( drds2 > 0. )
  {
    if( fabs(drds2 - 1.0) > 1.e-14 )
    {
      G4double normx = 1.0 / sqrt(drds2);
      for(G4int i=9;i<12;i++)  { vec[i] *= normx; }
    }
  }
}

