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
// $Id: G4MagIntegratorDriver.cc,v 1.57 2010-07-14 10:00:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// Implementation for class G4MagInt_Driver
// Tracking in space dependent magnetic field
//
// History of major changes:
//  8 Nov 01  J. Apostolakis:   Respect minimum step in AccurateAdvance
// 27 Jul 99  J. Apostolakis:   Ensured that AccurateAdvance does not loop 
//                              due to very small eps & step size (precision)
// 28 Jan 98  W. Wander:        Added ability for low order integrators
//  7 Oct 96  V. Grichine       First version
// --------------------------------------------------------------------

//#include "globals.hh"
//#include "G4GeometryTolerance.hh"
//#include <iomanip>

#include "GPMagIntegratorDriver.h"
#include "GPFieldTrack.h"
#include "GPUtils.h"

//  Stepsize can increase by no more than 5.0
//           and decrease by no more than 1/10. = 0.1
//
//const G4double G4MagInt_Driver::max_stepping_increase = 5.0;
//const G4double G4MagInt_Driver::max_stepping_decrease = 0.1;

//  The (default) maximum number of steps is Base
//  divided by the order of Stepper
//
//const G4int  G4MagInt_Driver::fMaxStepBase = 250;  // Was 5000

// ---------------------------------------------------------
//  Constructor
//
FQUALIFIER
void GPMagInt_Driver_Constructor(GPMagInt_Driver *This, 
				 G4double         hminimum, 
				 GPClassicalRK4  *pStepper,
				 G4int            numComponents,
				 G4int            statisticsVerbose)
{  
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12

  This->fSmallestFraction = 1.0e-12 ;
  This->fNoIntegrationVariables = numComponents; 
  This->fMinNoVars = 12;
  This->fNoVars =  GPimax(This->fNoIntegrationVariables, This->fMinNoVars);
  This->fStatisticsVerboseLevel = statisticsVerbose;
  This->fNoTotalSteps = 0; 
  This->fNoBadSteps = 0;
  This->fNoSmallSteps = 0;
  This->fNoInitialSmallSteps = 0; 
  This->fDyerr_max = 0.0;
  This->fDyerr_mx2 = 0.0;
  This->fDyerrPos_smTot = 0.0; 
  This->fDyerrPos_lgTot = 0.0;
  This->fDyerrVel_lgTot = 0.0; 
  This->fSumH_sm = 0.0;
  This->fSumH_lg = 0.0;
  This->fVerboseLevel = 0;

  GPMagInt_Driver_RenewStepperAndAdjust(This, pStepper );

  This->fMinimumStep = hminimum;
  This->fMaxNoSteps  = 
    This->fMaxStepBase/GPClassicalRK4_IntegratorOrder();

  if( (This->fVerboseLevel > 0) || (This->fStatisticsVerboseLevel > 1) )
  {
    ;
    //    G4cout << "MagIntDriver version: Accur-Adv: "
    //           << "invE_nS, QuickAdv-2sqrt with Statistics "
  }
}

// To add much printing for debugging purposes, uncomment the following
// and set verbose level to 1 or higher value !
// #define  G4DEBUG_FIELD 1    

// ---------------------------------------------------------

FQUALIFIER
G4bool GPMagInt_Driver_AccurateAdvance( GPMagInt_Driver *This,
					GPFieldTrack& y_current,
					G4double     hstep,
					G4double     eps,
					G4double hinitial )
{
  // Runge-Kutta driver with adaptive stepsize control. Integrate starting
  // values at y_current over hstep x2 with accuracy eps. 
  // On output ystart is replaced by values at the end of the integration 
  // interval. RightHandSide is the right-hand side of ODE system. 
  // The source is similar to odeint routine from NRC p.721-722 .

  G4int nstp, i, no_warnings=0;
  G4double x, hnext, hdid, h;

  G4double y[GPFieldTrack_ncompSVEC], dydx[GPFieldTrack_ncompSVEC];
  G4double ystart[GPFieldTrack_ncompSVEC], yEnd[GPFieldTrack_ncompSVEC]; 
  G4double  x1, x2;
  G4bool succeeded = true, lastStepSucceeded;

  G4double startCurveLength;

  G4int  noFullIntegr=0, noSmallIntegr = 0 ;
  //  static G4int  noGoodSteps =0 ;  // Bad = chord > curve-len 
  G4int  noGoodSteps =0 ;  // Bad = chord > curve-len 
  const  G4int  nvar= This->fNoVars;

  //  GPFieldTrack yStartFT(y_current);

  //  Ensure that hstep > 0
  //
  if( hstep <= 0.0 )
  { 
    if(hstep==0.0)
    {
      //      std::ostringstream message;
      //      message << "Proposed step is zero; hstep = " << hstep << " !";
      //      G4Exception("GPMagInt_Driver::AccurateAdvance()", 
      //                  "GeomField1001", JustWarning, message);
      return succeeded; 
    }
    else
    { 
      //      std::ostringstream message;
      //      message << "Invalid run condition." << G4endl
      //              << "Proposed step is negative; hstep = " << hstep << "." << G4endl
      //              << "Requested step cannot be negative! Aborting event.";
      //      G4Exception("GPMagInt_Driver::AccurateAdvance()", 
      //                  "GeomField0003", EventMustBeAborted, message);
      return false;
    }
  }

  GPFieldTrack_DumpToArray(&y_current, ystart );

  //  startCurveLength= y_current.GetCurveLength();
  startCurveLength= GPFieldTrack_GetCurveLength(&y_current);
  x1= startCurveLength; 
  x2= x1 + hstep;

  if ( (hinitial > 0.0) && (hinitial < hstep)
    && (hinitial > perMillion * hstep) )
  {
     h = hinitial;
  }
  else  //  Initial Step size "h" defaults to the full interval
  {
     h = hstep;
  }

  x = x1;

  for (i=0;i<nvar;i++)  { y[i] = ystart[i]; }

  G4bool lastStep= false;
  nstp=1;

  do
  {
    GPThreeVector StartPos = GPThreeVector_create( y[0], y[1], y[2] );

    // Old method - inline call to Equation of Motion
    //   pIntStepper->RightHandSide( y, dydx );
    // New method allows to cache field, or state (eg momentum magnitude)
    //    pIntStepper->ComputeRightHandSide( y, dydx );
    GPClassicalRK4_ComputeRightHandSide(This->pIntStepper, y, dydx );

    This->fNoTotalSteps++;

    // Perform the Integration
    //      
    if( h > This->fMinimumStep )
    { 
      GPMagInt_Driver_OneGoodStep(This, y,dydx,x,h,eps,hdid,hnext) ;
      //--------------------------------------
      lastStepSucceeded= (hdid == h);   
    }
    else
    {
      GPFieldTrack yFldTrk;
      GPFieldTrack_Constructor2(&yFldTrk, 
				GPThreeVector_create(0,0,0), 
				GPThreeVector_create(0,0,0), 
				0., 0., 0., 0.);

      G4double dchord_step, dyerr, dyerr_len;   // What to do with these ?
      //      yFldTrk.LoadFromArray(y, fNoIntegrationVariables); 
      //      yFldTrk.SetCurveLength( x );
      GPFieldTrack_LoadFromArray(&yFldTrk, y, This->fNoIntegrationVariables); 
      GPFieldTrack_SetCurveLength(&yFldTrk, x );

      GPMagInt_Driver_QuickAdvance(This, 
				   yFldTrk, dydx, h, dchord_step, dyerr_len ); 
      //-----------------------------------------------------

      GPFieldTrack_DumpToArray(&yFldTrk,y);    

      if( h == 0.0 )
      { 
	;
	//        G4Exception("GPMagInt_Driver::AccurateAdvance()",
	//                    "GeomField0003", FatalException,
	//                    "Integration Step became Zero!"); 
      }
      dyerr = dyerr_len / h;
      hdid= h;
      x += hdid;

      // Compute suggested new step
      hnext= GPMagInt_Driver_ComputeNewStepSize(This, dyerr/eps, h);

      // .. hnext= ComputeNewStepSize_WithinLimits( dyerr/eps, h);
      lastStepSucceeded= (dyerr<= eps);
    }

    if (lastStepSucceeded)  { noFullIntegr++; }
    else                    { noSmallIntegr++; }

    GPThreeVector EndPos = GPThreeVector_create( y[0], y[1], y[2] );

    // Check the endpoint
    G4double endPointDist= 
      GPThreeVector_mag(GPThreeVector_sub(EndPos,StartPos)); 
    if ( endPointDist >= hdid*(1.+perMillion) )
    {
      This->fNoBadSteps++;

      // Issue a warning only for gross differences -
      // we understand how small difference occur.
      if ( endPointDist >= hdid*(1.+perThousand) )
      { 
        no_warnings++;
      }
    }
    else
    {
      noGoodSteps ++;
    } 

    //  Avoid numerous small last steps
    if( (h < eps * hstep) || (h < This->fSmallestFraction * startCurveLength) )
    {
      // No more integration -- the next step will not happen
      lastStep = true;  
    }
    else
    {
      // Check the proposed next stepsize
      if(fabs(hnext) <= GPMagInt_Driver_Hmin(This))
      {
        // Make sure that the next step is at least Hmin.
        h = GPMagInt_Driver_Hmin(This);
      }
      else
      {
        h = hnext;
      }

      //  Ensure that the next step does not overshoot
      if ( x+h > x2 )
      {                // When stepsize overshoots, decrease it!
        h = x2 - x ;   // Must cope with difficult rounding-error
      }                // issues if hstep << x2

      if ( h == 0.0 )
      {
        // Cannot progress - accept this as last step - by default
        lastStep = true;
      }
    }
  } while ( ((nstp++)<=This->fMaxNoSteps) && (x < x2) && (!lastStep) );
     // Have we reached the end ?
     // --> a better test might be x-x2 > an_epsilon

  succeeded=  (x>=x2);  // If it was a "forced" last step

  for (i=0;i<nvar;i++)  { yEnd[i] = y[i]; }

  // Put back the values.
  //  y_current.LoadFromArray( yEnd, fNoIntegrationVariables );
  //  y_current.SetCurveLength( x );
  GPFieldTrack_LoadFromArray(&y_current, yEnd, This->fNoIntegrationVariables );
  GPFieldTrack_SetCurveLength(&y_current, x );

  if(nstp > This->fMaxNoSteps)
  {
    no_warnings++;
    succeeded = false;
  }


  return succeeded;
}  // end of AccurateAdvance ...........................

// ---------------------------------------------------------

FQUALIFIER
void GPMagInt_Driver_OneGoodStep( GPMagInt_Driver *This,
				  G4double y[],        // InOut
				  const G4double dydx[],
				  G4double& x,         // InOut
				  G4double htry,
				  G4double eps_rel_max,
				  G4double& hdid,      // Out
				  G4double& hnext )    // Out

// Driver for one Runge-Kutta Step with monitoring of local truncation error
// to ensure accuracy and adjust stepsize. Input are dependent variable
// array y[0,...,5] and its derivative dydx[0,...,5] at the
// starting value of the independent variable x . Also input are stepsize
// to be attempted htry, and the required accuracy eps. On output y and x
// are replaced by their new values, hdid is the stepsize that was actually
// accomplished, and hnext is the estimated next stepsize. 
// This is similar to the function rkqs from the book:
// Numerical Recipes in C: The Art of Scientific Computing (NRC), Second
// Edition, by William H. Press, Saul A. Teukolsky, William T.
// Vetterling, and Brian P. Flannery (Cambridge University Press 1992),
// 16.2 Adaptive StepSize Control for Runge-Kutta, p. 719

{
  G4double errmax_sq;
  G4double h, htemp, xnew ;

  G4double yerr[GPFieldTrack_ncompSVEC], ytemp[GPFieldTrack_ncompSVEC];

  h = htry ; // Set stepsize to the initial trial value

  G4double inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  G4double errpos_sq=0.0;    // square of displacement error
  G4double errvel_sq=0.0;    // square of momentum vector difference
  G4double errspin_sq=0.0;   // square of spin vector difference

  G4int iter;

  //  static G4int tot_no_trials=0; 
  G4int tot_no_trials=0; 
  const G4int max_trials=100; 

  GPThreeVector Spin = GPThreeVector_create(y[9],y[10],y[11]);
  G4bool     hasSpin= (GPThreeVector_mag2(Spin) > 0.0); 

  for (iter=0; iter<max_trials ;iter++)
  {
    tot_no_trials++;
    //@@@find who is drived from G4MagIntegratorStepper
    //    pIntStepper-> Stepper( y,dydx,h,ytemp,yerr); 
    //            *******
    GPClassicalRK4_Stepper(This->pIntStepper, 
			   y,dydx,h,ytemp,yerr); 

    G4double eps_pos = eps_rel_max * GPfmax(h, This->fMinimumStep); 
    G4double inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos); 

    // Evaluate accuracy
    //
    errpos_sq =  GPsqr(yerr[0]) + GPsqr(yerr[1]) + GPsqr(yerr[2]) ;
    errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

    // Accuracy for momentum
    errvel_sq =  (GPsqr(yerr[3]) + GPsqr(yerr[4]) + GPsqr(yerr[5]) )
               / (GPsqr(y[3]) + GPsqr(y[4]) + GPsqr(y[5]) );
    errvel_sq *= inv_eps_vel_sq;
    errmax_sq = GPfmax( errpos_sq, errvel_sq ); // Square of maximum error

    if( hasSpin )
    { 
      // Accuracy for spin
      errspin_sq =  ( GPsqr(yerr[9]) + GPsqr(yerr[10]) + GPsqr(yerr[11]) )
                 /  ( GPsqr(y[9]) + GPsqr(y[10]) + GPsqr(y[11]) );
      errspin_sq *= inv_eps_vel_sq;
      errmax_sq = GPfmax( errmax_sq, errspin_sq ); 
   }

    if ( errmax_sq <= 1.0 )  { break; } // Step succeeded. 

    // Step failed; compute the size of retrial Step.
    htemp = GPMagInt_Driver_GetSafety(This)*h* 
      pow( errmax_sq, 0.5*GPMagInt_Driver_GetPshrnk(This) );

    if (htemp >= 0.1*h)  { h = htemp; }  // Truncation error too large,
    else  { h = 0.1*h; }                 // reduce stepsize, but no more
                                         // than a factor of 10
    xnew = x + h;
    if(xnew == x)
    {
      //      G4cerr << "G4MagIntegratorDriver::OneGoodStep:" << G4endl
      //             << "  Stepsize underflow in Stepper " << G4endl ;
      //      G4cerr << "  Step's start x=" << x << " and end x= " << xnew 
      //             << " are equal !! " << G4endl
      //             <<"  Due to step-size= " << h 
      //             << " . Note that input step was " << htry << G4endl;
      break;
    }
  }

  // Compute size of next Step
  if (errmax_sq > This->errcon*This->errcon)
  { 
    hnext = GPMagInt_Driver_GetSafety(This)*h*
      pow(errmax_sq, 0.5*GPMagInt_Driver_GetPgrow(This));
  }
  else
  {
    hnext = max_stepping_increase*h ; // No more than a factor of 5 increase
  }
  x += (hdid = h);

  for(G4int k=0 ; k < This->fNoIntegrationVariables ; k++) { y[k] = ytemp[k]; }

  return;
}   // end of  OneGoodStep .............................

//----------------------------------------------------------------------

// QuickAdvance just tries one Step - it does not ensure accuracy
//
FQUALIFIER
G4bool  GPMagInt_Driver_QuickAdvance2(GPMagInt_Driver *This,       
				      GPFieldTrack& y_posvel,         // INOUT
				      const G4double     dydx[],  
				      G4double     hstep,       // In
				      G4double&    dchord_step,
				      G4double&    dyerr_pos_sq,
				      G4double&    dyerr_mom_rel_sq )  
{
  //  G4Exception("GPMagInt_Driver::QuickAdvance()", "GeomField0001",
  //              FatalException, "Not yet implemented."); 

  // Use the parameters of this method, to please compiler
  dchord_step = dyerr_pos_sq = hstep * hstep * dydx[0]; 
  dyerr_mom_rel_sq = GPThreeVector_mag2(GPFieldTrack_GetPosition(&y_posvel));
  return true;
}

//----------------------------------------------------------------------
FQUALIFIER
G4bool  GPMagInt_Driver_QuickAdvance(GPMagInt_Driver *This,        
				     GPFieldTrack& y_posvel,         // INOUT
				     const G4double     dydx[],  
				     G4double     hstep,       // In
				     G4double&    dchord_step,
				     G4double&    dyerr )
{
  G4double dyerr_pos_sq, dyerr_mom_rel_sq;  
  G4double yerr_vec[GPFieldTrack_ncompSVEC],
           yarrin[GPFieldTrack_ncompSVEC], yarrout[GPFieldTrack_ncompSVEC]; 
  G4double s_start;
  G4double dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

  //  static 
  G4int no_call=0; 
  no_call ++; 

  // Move data into array
  GPFieldTrack_DumpToArray(&y_posvel, yarrin );      //  yarrin  <== y_posvel 
  s_start = GPFieldTrack_GetCurveLength(&y_posvel);

  // Do an Integration Step
  //  pIntStepper-> Stepper(yarrin, dydx, hstep, yarrout, yerr_vec) ; 
  GPClassicalRK4_Stepper(This->pIntStepper, 
			 yarrin, dydx, hstep, yarrout, yerr_vec) ; 
  //            *******

  // Estimate curve-chord distance
  //  dchord_step= G4MagIntegratorStepper_DistChord(pIntStepper);
  dchord_step= GPClassicalRK4_DistChord(This->pIntStepper);
  //                         *********

  // Put back the values.  yarrout ==> y_posvel
  GPFieldTrack_LoadFromArray(&y_posvel, yarrout, This->fNoIntegrationVariables);
  GPFieldTrack_SetCurveLength(&y_posvel, s_start + hstep );

  // A single measure of the error   
  //      TO-DO :  account for  energy,  spin, ... ? 
  vel_mag_sq   = ( GPsqr(yarrout[3])+GPsqr(yarrout[4])+GPsqr(yarrout[5]) );
  inv_vel_mag_sq = 1.0 / vel_mag_sq; 
  dyerr_pos_sq = ( GPsqr(yerr_vec[0])+GPsqr(yerr_vec[1])+GPsqr(yerr_vec[2]));
  dyerr_mom_sq = ( GPsqr(yerr_vec[3])+GPsqr(yerr_vec[4])+GPsqr(yerr_vec[5]));
  dyerr_mom_rel_sq = dyerr_mom_sq * inv_vel_mag_sq;

  // Calculate also the change in the momentum squared also ???
  // G4double veloc_square = y_posvel.GetVelocity().mag2();
  // ...

  if( dyerr_pos_sq > ( dyerr_mom_rel_sq * GPsqr(hstep) ) )
  {
    dyerr = sqrt(dyerr_pos_sq);
  }
  else
  {
    // Scale it to the current step size - for now
    dyerr = sqrt(dyerr_mom_rel_sq) * hstep;
  }

  return true;
}

// --------------------------------------------------------------------------

//  This method computes new step sizes - but does not limit changes to
//   within  certain factors
// 
FQUALIFIER
G4double GPMagInt_Driver_ComputeNewStepSize(GPMagInt_Driver *This, 
                          G4double  errMaxNorm,    // max error  (normalised)
                          G4double  hstepCurrent)  // current step size
{
  G4double hnew;

  // Compute size of next Step for a failed step
  if(errMaxNorm > 1.0 )
  {
    // Step failed; compute the size of retrial Step.
    hnew = GPMagInt_Driver_GetSafety(This)*hstepCurrent*
      pow(errMaxNorm,GPMagInt_Driver_GetPshrnk(This)) ;
  } else if(errMaxNorm > 0.0 ) {
    // Compute size of next Step for a successful step
    hnew = GPMagInt_Driver_GetSafety(This)*hstepCurrent*
      pow(errMaxNorm,GPMagInt_Driver_GetPgrow(This)) ;
  } else {
    // if error estimate is zero (possible) or negative (dubious)
    hnew = max_stepping_increase * hstepCurrent; 
  }

  return hnew;
}

// ---------------------------------------------------------------------------

// This method computes new step sizes limiting changes within certain factors
// 
// It shares its logic with AccurateAdvance.
// They are kept separate currently for optimisation.
//
FQUALIFIER
G4double GPMagInt_Driver_ComputeNewStepSize_WithinLimits(GPMagInt_Driver *This, 
			  G4double  errMaxNorm,    // max error  (normalised)
                          G4double  hstepCurrent)  // current step size
{
  G4double hnew;

  // Compute size of next Step for a failed step
  if (errMaxNorm > 1.0 )
  {
    // Step failed; compute the size of retrial Step.
    hnew = GPMagInt_Driver_GetSafety(This)*hstepCurrent*
      pow(errMaxNorm,GPMagInt_Driver_GetPshrnk(This)) ;
  
    if (hnew < max_stepping_decrease*hstepCurrent)
    {
      hnew = max_stepping_decrease*hstepCurrent ;
                         // reduce stepsize, but no more
                         // than this factor (value= 1/10)
    }
  }
  else
  {
    // Compute size of next Step for a successful step
    if (errMaxNorm > This->errcon)
     { hnew = GPMagInt_Driver_GetSafety(This)*hstepCurrent*
	 pow(errMaxNorm,GPMagInt_Driver_GetPgrow(This)); }
    else  // No more than a factor of 5 increase
     { hnew = max_stepping_increase * hstepCurrent; }
  }
  return hnew;
}

// ---------------------------------------------------------------------------

FQUALIFIER
void GPMagInt_Driver_SetSmallestFraction(GPMagInt_Driver *This,
					 G4double newFraction)
{
  if( (newFraction > 1.e-16) && (newFraction < 1e-8) )
  {
    This->fSmallestFraction= newFraction;
  }
  else
  { 
    //    G4cerr << "Warning: SmallestFraction not changed. " << G4endl
    //           << "  Proposed value was " << newFraction << G4endl
    //           << "  Value must be between 1.e-8 and 1.e-16" << G4endl;
  }
}

//inline functions
//G4MagIntegratorDriver.icc

FQUALIFIER
G4double GPMagInt_Driver_GetHmin( GPMagInt_Driver *This )
{
  return This->fMinimumStep;
} 

FQUALIFIER
G4double GPMagInt_Driver_Hmin( GPMagInt_Driver *This )
{
  return This->fMinimumStep;
}

FQUALIFIER
G4double GPMagInt_Driver_GetSafety( GPMagInt_Driver *This )
{
  return This->safety;
}

FQUALIFIER
G4double GPMagInt_Driver_GetPshrnk( GPMagInt_Driver *This )
{
  return This->pshrnk;
} 

FQUALIFIER
G4double GPMagInt_Driver_GetPgrow( GPMagInt_Driver *This )
{
  return This->pgrow;
}
 
FQUALIFIER
G4double GPMagInt_Driver_GetErrcon( GPMagInt_Driver *This )
{
  return This->errcon;
}

FQUALIFIER
void GPMagInt_Driver_SetHmin( GPMagInt_Driver *This, G4double newval)
{
  This->fMinimumStep = newval;
} 
FQUALIFIER
G4double GPMagInt_Driver_ComputeAndSetErrcon( GPMagInt_Driver *This )
{
  This->errcon = pow(max_stepping_increase/GPMagInt_Driver_GetSafety(This),
		    1.0/GPMagInt_Driver_GetPgrow(This));
  return This->errcon;
} 

FQUALIFIER
void GPMagInt_Driver_ReSetParameters( GPMagInt_Driver *This,
				      G4double new_safety)
{
  This->safety = new_safety;
  This->pshrnk = -1.0 /GPClassicalRK4_IntegratorOrder();
  This->pgrow  = -1.0 /(1.0+GPClassicalRK4_IntegratorOrder());
  GPMagInt_Driver_ComputeAndSetErrcon(This);
}

FQUALIFIER
void GPMagInt_Driver_SetSafety( GPMagInt_Driver *This, G4double val)
{ 
  This->safety=val;
  GPMagInt_Driver_ComputeAndSetErrcon(This);
}

FQUALIFIER
void GPMagInt_Driver_SetPgrow( GPMagInt_Driver *This, G4double  val)
{ 
  This->pgrow=val;
  GPMagInt_Driver_ComputeAndSetErrcon(This); 
}

FQUALIFIER
void GPMagInt_Driver_SetErrcon( GPMagInt_Driver *This, G4double val)
{ 
  This->errcon=val;
}

FQUALIFIER
void GPMagInt_Driver_RenewStepperAndAdjust( GPMagInt_Driver *This,
					    GPClassicalRK4 *pItsStepper)
{  
  This->pIntStepper = pItsStepper; 
  GPMagInt_Driver_ReSetParameters(This,0.0);
}
FQUALIFIER
void GPMagInt_Driver_SetChargeMomentumMass( GPMagInt_Driver *This,
					    G4double particleCharge, // e+
					    G4double MomentumXc,
					    G4double Mass )
{ 
  GPEquationOfMotion_SetChargeMomentumMass(
		     GPClassicalRK4_GetEquationOfMotion(This->pIntStepper),
		     particleCharge, MomentumXc, Mass); 
}

FQUALIFIER
const GPClassicalRK4* GPMagInt_Driver_GetStepper(GPMagInt_Driver *This)
{
  return This->pIntStepper;
}

FQUALIFIER
G4int GPMagInt_Driver_GetMaxNoSteps( GPMagInt_Driver *This )
{
  return This->fMaxNoSteps;
}

FQUALIFIER
void GPMagInt_Driver_SetMaxNoSteps( GPMagInt_Driver *This, G4int val)
{
  This->fMaxNoSteps= val;
}

FQUALIFIER
void GPMagInt_Driver_GetDerivatives( GPMagInt_Driver *This,
				     GPFieldTrack &y_curr,// const, INput
				     G4double     dydx[])  // OUTput
{ 
  G4double  tmpValArr[GPFieldTrack_ncompSVEC];
  GPFieldTrack_DumpToArray(&y_curr, tmpValArr  );
  GPClassicalRK4_RightHandSide(This->pIntStepper, tmpValArr , dydx );
}

FQUALIFIER
G4double GPMagInt_Driver_GetVerboseLevel( GPMagInt_Driver *This )
{
  return This->fVerboseLevel;
} 

FQUALIFIER 
void GPMagInt_Driver_SetVerboseLevel( GPMagInt_Driver *This, G4int newLevel)
{
  This->fVerboseLevel= newLevel;
}

FQUALIFIER
G4double GPMagInt_Driver_GetSmallestFraction( GPMagInt_Driver *This )
{
  return This->fSmallestFraction; 
} 
