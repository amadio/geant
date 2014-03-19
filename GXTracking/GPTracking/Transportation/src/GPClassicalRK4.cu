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
// $Id: G4MagIntegratorStepper.cc,v 1.12 2009-11-05 18:31:15 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------

#include "GPClassicalRK4.h"
#include "GPUtils.h"
#include "GPLineSection.h"

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
void GPClassicalRK4_Constructor( GPClassicalRK4 *This, 
				 GPEquationOfMotion* EqRhs,
				 G4int numberOfVariables)
{

  //  This->fNoIntegrationVariables = numberOfVariables;
  //  This->fNoStateVariables = numberOfVariables;

  GPClassicalRK4_G4MagErrorStepper_Constructor(This, 
  					       EqRhs,numberOfVariables,0);

  int noVariables= GPimax(numberOfVariables,8); // For Time .. 7+1

  //  dydxm = new G4double[noVariables];
  //  dydxt = new G4double[noVariables]; 
  //  yt    = new G4double[noVariables]; 

  //only cuda compatible - not good for both host and device! 
  //  cudaMalloc((void**)&(This->dydxm), noVariables*sizeof(G4double));
  //  cudaMalloc((void**)&(This->dydxt), noVariables*sizeof(G4double));
  //  cudaMalloc((void**)&(This->yt),    noVariables*sizeof(G4double));

}

// G4ClassicalRK4::~G4ClassicalRK4()
FQUALIFIER
void GPClassicalRK4_Destructor( GPClassicalRK4 *This) 
{ 
  //only cuda compatible
  //  cudaFree(This->dydxm);
  //  cudaFree(This->dydxt);
  //  cudaFree(This->yt);

  GPClassicalRK4_G4MagErrorStepper_Destructor(This);

}

// G4ClassicalRK4::IntegratorOrder
FQUALIFIER
G4int GPClassicalRK4_IntegratorOrder() { return 4; }

// G4ClassicalRK4::DumbStepper
FQUALIFIER
void GPClassicalRK4_DumbStepper( GPClassicalRK4 *This,
				 const G4double  yIn[],
				 const G4double  dydx[],
				 G4double  h,
				 G4double  yOut[])
{
  const G4int nvar = 
    GPClassicalRK4_GetNumberOfVariables(This); //fNumberOfVariables(); 

  G4int i;
  G4double  hh = h*0.5 , h6 = h/6.0  ;

  // Initialise time to t0, needed when it is not updated by the integration.
  //        [ Note: Only for time dependent fields (usually electric) 
  //                  is it neccessary to integrate the time.] 
  This->yt[7]   = yIn[7]; 
  yOut[7] = yIn[7];

  // 1st Step K1=h*dydx
  for(i=0;i<nvar;i++)
  {
    This->yt[i] = yIn[i] + hh*dydx[i] ;     
  }

  // 2nd Step K2=h*dydxt
  GPClassicalRK4_RightHandSide(This, This->yt,This->dydxt) ; 

  for(i=0;i<nvar;i++)
  { 
    This->yt[i] = yIn[i] + hh*This->dydxt[i] ;
  }

  // 3rd Step K3=h*dydxm
  GPClassicalRK4_RightHandSide(This, This->yt,This->dydxm) ;  

  // now dydxm=(K2+K3)/h
  for(i=0;i<nvar;i++)
  {
    This->yt[i]   = yIn[i] + h*This->dydxm[i] ;
    This->dydxm[i] += This->dydxt[i] ;                 
  }

  // 4th Step K4=h*dydxt
  GPClassicalRK4_RightHandSide(This, This->yt,This->dydxt) ;  
 
  // Final RK4 output
  for(i=0;i<nvar;i++)    
  {
    //+K1/6+K4/6+(K2+K3)/3
    yOut[i] = yIn[i]+h6*(dydx[i]+This->dydxt[i]+2.0*This->dydxm[i]); 
  }
  if ( nvar == 12 )  
  { 
    GPClassicalRK4_NormalisePolarizationVector (This, yOut ); 
  }
  
}  // end of DumbStepper ....................................................

// G4ClassicalRK4::StepWithEst
/*
FQUALIFIER
void GPClassicalRK4_StepWithEst( const G4double*,
				 const G4double*,
				 G4double,
				 G4double*,
				 G4double&,
				 G4double&,
				 const G4double*,
				 G4double*  ) 
{
  ;
  //  G4Exception("G4ClassicalRK4::StepWithEst()", "GeomField0001",
  //              FatalException, "Method no longer used.");

} 
*/

//--------------------------------------------------------
// class G4MagErrorStepper : public G4MagIntegratorStepper
//--------------------------------------------------------  

// G4MagErrorStepper::G4MagErrorStepper
FQUALIFIER
void GPClassicalRK4_G4MagErrorStepper_Constructor(GPClassicalRK4 *This,
				     GPEquationOfMotion *EquationRhs,
                                     G4int numberOfVariables, 
                                     G4int numStateVariables)
{

  GPClassicalRK4_G4MagIntegratorStepper_Constructor(This,
		       EquationRhs,numberOfVariables,numStateVariables);

  G4int nvar = GPimax(GPClassicalRK4_GetNumberOfVariables(This), 8);
  //  yMiddle=     new G4double[nvar]; 
  //  dydxMid=     new G4double[nvar];
  //  yInitial=    new G4double[nvar];
  //  yOneStep=    new G4double[nvar];

  //only cuda compatible
  //  cudaMalloc((void**)&(This->yMiddle),  nvar*sizeof(G4double));
  //  cudaMalloc((void**)&(This->dydxMid),  nvar*sizeof(G4double));
  //  cudaMalloc((void**)&(This->yInitial), nvar*sizeof(G4double));
  //  cudaMalloc((void**)&(This->yOneStep), nvar*sizeof(G4double));

}

// G4MagErrorStepper::~G4MagErrorStepper()
FQUALIFIER
void GPClassicalRK4_G4MagErrorStepper_Destructor(GPClassicalRK4 *This)
{
  //   delete[] yMiddle;
  //   delete[] dydxMid;
  //   delete[] yInitial;
  //   delete[] yOneStep;

  //  cudaFree(This->yMiddle);
  //  cudaFree(This->dydxMid);
  //  cudaFree(This->yInitial);
  //  cudaFree(This->yOneStep);
}

// G4MagErrorStepper::Stepper
FQUALIFIER
void GPClassicalRK4_Stepper( GPClassicalRK4 *This, 
			     const G4double yInput[],
			     const G4double dydx[],
			     G4double hstep,
			     G4double yOutput[],
			     G4double yError []      )
{  
  const G4int nvar = GPClassicalRK4_GetNumberOfVariables(This) ;
  const G4int maxvar= GPClassicalRK4_GetNumberOfStateVariables(This);

  G4int i;
  // correction for Richardson Extrapolation.
  G4double correction = 1./((1 << GPClassicalRK4_IntegratorOrder()) -1 );
  
  //  Saving yInput because yInput and yOutput can be aliases for same array
  
  for(i=0;i<nvar;i++) This->yInitial[i]=yInput[i];
  This->yInitial[7]= yInput[7];//Copy the time in case even if not really needed
  This->yMiddle[7] = yInput[7]; //Copy the time from initial value 
  This->yOneStep[7] = yInput[7];//As it contributes to final value of yOutput ?
  // yOutput[7] = yInput[7];  // -> dumb stepper does it too for RK4
  for(i=nvar;i<maxvar;i++) yOutput[i]=yInput[i];
  // yError[7] = 0.0;         
  
  G4double halfStep = hstep * 0.5; 
  
  // Do two half steps
  
  GPClassicalRK4_DumbStepper(This, This->yInitial, dydx,
			     halfStep, This->yMiddle);
  GPClassicalRK4_RightHandSide(This, This->yMiddle, This->dydxMid);    
  GPClassicalRK4_DumbStepper(This, This->yMiddle, This->dydxMid, 
			     halfStep, yOutput); 
  
  // Store midpoint, chord calculation
  This->fMidPoint = GPThreeVector_create( This->yMiddle[0],  
					  This->yMiddle[1],  
					  This->yMiddle[2]); 
  
  // Do a full Step
  GPClassicalRK4_DumbStepper(This, This->yInitial, dydx, hstep, This->yOneStep);
  for(i=0;i<nvar;i++) {
    yError [i] = yOutput[i] - This->yOneStep[i] ;
    yOutput[i] += yError[i]*correction ;  // Provides accuracy increased
    // by 1 order via the 
    // Richardson Extrapolation  
  }
  
  This->fInitialPoint = GPThreeVector_create( This->yInitial[0], 
					      This->yInitial[1], 
					      This->yInitial[2]); 

  This->fFinalPoint   = GPThreeVector_create( yOutput[0],  
					      yOutput[1],  
					      yOutput[2]); 
  
  return ;
}

// G4MagErrorStepper::DistChord
FQUALIFIER
G4double GPClassicalRK4_DistChord(GPClassicalRK4 *This)
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

  if (GPThreeVector_nequal(This->fInitialPoint,This->fFinalPoint)) {
     distLine= GPLineSection_Distline( This->fMidPoint, This->fInitialPoint, 
				       This->fFinalPoint );
     // This is a class method that gives distance of Mid 
     //  from the Chord between the Initial and Final points.

     distChord = distLine;
  }else{
    distChord = GPThreeVector_mag(GPThreeVector_sub(This->fMidPoint,
						    This->fInitialPoint));
  }

  return distChord;
}

//--------------------------------------------------------------
// class G4MagIntegratorStepper
//--------------------------------------------------------------  

// G4MagIntegratorStepper::G4MagIntegratorStepper
FQUALIFIER
void GPClassicalRK4_G4MagIntegratorStepper_Constructor(GPClassicalRK4 *This,
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
void GPClassicalRK4_ComputeRightHandSide( GPClassicalRK4 *This,
					  const G4double y[], 
					  G4double dydx[] ) 
{
  GPClassicalRK4_RightHandSide(This, y, dydx );
}

// G4MagIntegratorStepper::GetEquationOfMotion
FQUALIFIER
GPEquationOfMotion* GPClassicalRK4_GetEquationOfMotion(GPClassicalRK4 *This)
{
  return  This->fEquation_Rhs;
} 

// G4MagIntegratorStepper::SetEquationOfMotion
FQUALIFIER
void GPClassicalRK4_SetEquationOfMotion(GPClassicalRK4 *This,
					GPEquationOfMotion* newEquation)
{
  if( newEquation != 0 )
  {
    This->fEquation_Rhs= newEquation;
  }
} 

// G4MagIntegratorStepper::GetNumberOfVariables
FQUALIFIER
G4int GPClassicalRK4_GetNumberOfVariables(GPClassicalRK4 *This) 
{
  return This->fNoIntegrationVariables;
}

// G4MagIntegratorStepper::GetNumberOfStateVariables
FQUALIFIER
G4int GPClassicalRK4_GetNumberOfStateVariables(GPClassicalRK4 *This) 
{
  return This->fNoStateVariables;
}

// G4MagIntegratorStepper::RightHandSide
FQUALIFIER 
void GPClassicalRK4_RightHandSide(GPClassicalRK4 *This,
				  const  G4double y[], 
				  G4double dydx[] )   
{
  GPEquationOfMotion_RightHandSide(This->fEquation_Rhs, y, dydx);
}

// G4MagIntegratorStepper::NormaliseTangentVector
FQUALIFIER 
void GPClassicalRK4_NormaliseTangentVector(GPClassicalRK4 *This, 
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
void GPClassicalRK4_NormalisePolarizationVector(GPClassicalRK4 *This, 
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

