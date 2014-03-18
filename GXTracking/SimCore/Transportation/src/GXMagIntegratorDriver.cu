#include "GXMagIntegratorDriver.h"
#include "GXFieldTrack.h"
#include "GPUtils.h"

#include "stdio.h"

FQUALIFIER
void GXMagInt_Driver_Constructor(GXMagInt_Driver *This, 
				 G4double         hminimum, 
				 GXClassicalRK4  *pStepper)
{  
  // assume that the integration order = 4 and the number of variables = 6

  This->fMinimumStep = hminimum;
  This->pIntStepper = pStepper; 

  This->safety =  0.9;
  This->pshrnk = -0.25; // -1.0/(integration order = 4)
  This->pgrow  = -0.20; // -1.0/(integration order + 1)

  This->errcon = pow(max_stepping_increase/This->safety, 1.0/This->pgrow);

  This->fSmallestFraction = 1.0e-12 ;

  This->fMaxNoSteps = 62 ; //(fMaxStepBase=250)/(integration order=4);

}

// ---------------------------------------------------------

FQUALIFIER
G4bool GXMagInt_Driver_AccurateAdvance( GXMagInt_Driver *This,
					GXFieldTrack& y_current,
					G4double     hstep,
					G4double     eps,
					G4double hinitial )
{
  // Runge-Kutta driver with adaptive stepsize control. Integrate starting
  // values at y_current over hstep x2 with accuracy eps. 
  // On output ystart is replaced by values at the end of the integration 
  // interval. RightHandSide is the right-hand side of ODE system. 
  // The source is similar to odeint routine from NRC p.721-722 .

  G4int nstp, i =0;
  G4double x, hnext, hdid, h;

  G4double y[GXFieldTrack_ncompSVEC], dydx[GXFieldTrack_ncompSVEC];
  G4double ystart[GXFieldTrack_ncompSVEC], yEnd[GXFieldTrack_ncompSVEC]; 
  G4double  x1, x2;
  G4bool succeeded = true;

  G4double startCurveLength;

  //  Ensure that hstep > 0
  //
  if( hstep <= 0.0 )
  { 
    if(hstep==0.0)
    {
      //warning message
      return succeeded; 
    }
    else
    { 
      //exception
      return false;
    }
  }

  GXFieldTrack_DumpToArray(&y_current, ystart );

  //  startCurveLength= y_current.GetCurveLength();
  startCurveLength= GXFieldTrack_GetCurveLength(&y_current);
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

  for (i=0 ;i < GXFieldTrack_ncompSVEC ;i++)  { y[i] = ystart[i]; }

  G4bool lastStep= false;
  nstp=1;

  do
  {
    GPThreeVector StartPos = GPThreeVector_create( y[0], y[1], y[2] );

    GXClassicalRK4_RightHandSide(This->pIntStepper, y, dydx );

    // Perform the Integration
    //      
    if( h > This->fMinimumStep )
    { 
      GXMagInt_Driver_OneGoodStep(This, y,dydx,x,h,eps,hdid,hnext) ;
    }
    else
    {
      GXFieldTrack yFldTrk;
      GXFieldTrack_Constructor(&yFldTrk, 
			       GPThreeVector_create(0,0,0), 
			       GPThreeVector_create(0,0,0), 
			       0., 0., 0., 0.);

      G4double dchord_step, dyerr, dyerr_len;   // What to do with these ?
      GXFieldTrack_LoadFromArray(&yFldTrk, y); 
      GXFieldTrack_SetCurveLength(&yFldTrk, x );

      GXMagInt_Driver_QuickAdvance(This, 
				   yFldTrk, dydx, h, dchord_step, dyerr_len ); 
      //-----------------------------------------------------

      GXFieldTrack_DumpToArray(&yFldTrk,y);    

      if( h == 0.0 )
      { 
	; //invoke exception
      }

      dyerr = dyerr_len / h;
      hdid= h;
      x += hdid;

      // Compute suggested new step
      hnext= GXMagInt_Driver_ComputeNewStepSize(This, dyerr/eps, h);
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
      if(fabs(hnext) <= GXMagInt_Driver_GetHmin(This))
      {
        // Make sure that the next step is at least Hmin.
        h = GXMagInt_Driver_GetHmin(This);
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

  for (i=0; i < GXFieldTrack_ncompSVEC ;i++)  { yEnd[i] = y[i]; }

  // Put back the values.
  GXFieldTrack_LoadFromArray(&y_current, yEnd );
  GXFieldTrack_SetCurveLength(&y_current, x );

  if(nstp > This->fMaxNoSteps)
  {
    succeeded = false;
  }

  return succeeded;
}  // end of AccurateAdvance ...........................

// ---------------------------------------------------------

FQUALIFIER
void GXMagInt_Driver_OneGoodStep( GXMagInt_Driver *This,
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

  G4double yerr[GXFieldTrack_ncompSVEC], ytemp[GXFieldTrack_ncompSVEC];

  h = htry ; // Set stepsize to the initial trial value

  G4double inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  G4double errpos_sq=0.0;    // square of displacement error
  G4double errvel_sq=0.0;    // square of momentum vector difference

  const G4int max_trials=100; 

  G4int ntest = 0;
  for (G4int iter=0; iter<max_trials ;iter++)
  {
    ntest++;
    GXClassicalRK4_Stepper(This->pIntStepper,y,dydx,h,ytemp,yerr); 

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

    //    printf("errpos_sq errvel_sq errmax_sq = %f %f %f\n",
    //	   errpos_sq,errvel_sq,errmax_sq);

    if ( errmax_sq <= 1.0 )  { break; } // Step succeeded. 

    // Step failed; compute the size of retrial Step.
    htemp = GXMagInt_Driver_GetSafety(This)*h* 
      pow( errmax_sq, 0.5*GXMagInt_Driver_GetPshrnk(This) );

    if (htemp >= 0.1*h)  { h = htemp; }  // Truncation error too large,
    else  { h = 0.1*h; }                 // reduce stepsize, but no more
                                         // than a factor of 10
    xnew = x + h;
    if(xnew == x)
    {
      //error message
      break;
    }
  }

  // Compute size of next Step
  if (errmax_sq > This->errcon*This->errcon)
  { 
    hnext = GXMagInt_Driver_GetSafety(This)*h*
      pow(errmax_sq, 0.5*GXMagInt_Driver_GetPgrow(This));
  }
  else
  {
    hnext = max_stepping_increase*h ; // No more than a factor of 5 increase
  }
  x += (hdid = h);

  for(G4int k=0 ; k < GXFieldTrack_ncompSVEC ; k++) { y[k] = ytemp[k]; }

  //  printf("number of loops = %d\n",ntest);

}   // end of  OneGoodStep .............................


//@@@G4FWP test implemenation and should be used only for performance evaluation
FQUALIFIER
void GXMagInt_Driver_OneGoodStep2( GXMagInt_Driver *This,
				   G4double y[],        // InOut
				   const G4double dydx[],
				   G4double& x,         // InOut
				   G4double htry,
				   G4double eps_rel_max,
				   G4double& hdid,      // Out
				   G4double& hnext )    // Out
{
  //  G4double errmax_sq;
  //  G4double h, htemp, xnew ;
  G4double h;

  G4double yerr[GXFieldTrack_ncompSVEC], ytemp[GXFieldTrack_ncompSVEC];

  h = htry ; // Set stepsize to the initial trial value

  G4double inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  G4double errpos_sq=0.0;    // square of displacement error
  G4double errvel_sq=0.0;    // square of momentum vector difference

  //  const G4int max_trials=100; 

  //  G4int ntest = 0;
  //  for (G4int iter=0; iter<max_trials ;iter++)
  {
    //    ntest++;
    GXClassicalRK4_Stepper(This->pIntStepper,y,dydx,h,ytemp,yerr); 

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
    //    errmax_sq = GPfmax( errpos_sq, errvel_sq ); // Square of maximum error

    //    printf("errpos_sq errvel_sq errmax_sq = %f %f %f\n",
    //	   errpos_sq,errvel_sq,errmax_sq);

    //    if ( errmax_sq <= 1.0 )  { break; } // Step succeeded. 

    // Step failed; compute the size of retrial Step.
    //    htemp = GXMagInt_Driver_GetSafety(This)*h* 
    //      pow( errmax_sq, 0.5*GXMagInt_Driver_GetPshrnk(This) );

    //    if (htemp >= 0.1*h)  { h = htemp; }  // Truncation error too large,
    //    else  { h = 0.1*h; }                 // reduce stepsize, but no more
                                         // than a factor of 10
    //    xnew = x + h;
    //    if(xnew == x)
    //    {
      //error message
    //      break;
    //    }
  }

  // Compute size of next Step
  /* 
  if (errmax_sq > This->errcon*This->errcon)
  { 
    hnext = GXMagInt_Driver_GetSafety(This)*h*
      pow(errmax_sq, 0.5*GXMagInt_Driver_GetPgrow(This));
  }
  else
  {
    hnext = max_stepping_increase*h ; // No more than a factor of 5 increase
  }
  x += (hdid = h);
  */

  for(G4int k=0 ; k < GXFieldTrack_ncompSVEC ; k++) { y[k] = ytemp[k]; }

  //  printf("number of loops = %d\n",ntest);

}   // end of  OneGoodStep .............................


//----------------------------------------------------------------------
FQUALIFIER
G4bool  GXMagInt_Driver_QuickAdvance(GXMagInt_Driver *This,        
				     GXFieldTrack& y_posvel,         // INOUT
				     const G4double     dydx[],  
				     G4double     hstep,       // In
				     G4double&    dchord_step,
				     G4double&    dyerr )
{
  G4double dyerr_pos_sq, dyerr_mom_rel_sq;  
  G4double yerr_vec[GXFieldTrack_ncompSVEC],
           yarrin[GXFieldTrack_ncompSVEC], yarrout[GXFieldTrack_ncompSVEC]; 
  G4double s_start;
  G4double dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

  // Move data into array
  GXFieldTrack_DumpToArray(&y_posvel, yarrin );      //  yarrin  <== y_posvel 
  s_start = GXFieldTrack_GetCurveLength(&y_posvel);

  // Do an Integration Step
  GXClassicalRK4_Stepper(This->pIntStepper, 
			 yarrin, dydx, hstep, yarrout, yerr_vec) ; 


   // Estimate curve-chord distance
  dchord_step= GXClassicalRK4_DistChord(This->pIntStepper);
  //                         *********

  // Put back the values.  yarrout ==> y_posvel
  GXFieldTrack_LoadFromArray(&y_posvel, yarrout);
  GXFieldTrack_SetCurveLength(&y_posvel, s_start + hstep );

  // A single measure of the error   
  //      TO-DO :  account for  energy,  spin, ... ? 


  vel_mag_sq = (yarrout[3]*yarrout[3] + yarrout[4]*yarrout[4] +
		yarrout[5]*yarrout[5] );
  inv_vel_mag_sq = 1.0 / vel_mag_sq; 
  dyerr_pos_sq = ( yerr_vec[0]*yerr_vec[0]+yerr_vec[1]*yerr_vec[1]+
		   yerr_vec[2]*yerr_vec[2]);
  dyerr_mom_sq = ( yerr_vec[3]*yerr_vec[3]+yerr_vec[4]*yerr_vec[4]+
		   yerr_vec[5]*yerr_vec[5]);
  dyerr_mom_rel_sq = dyerr_mom_sq * inv_vel_mag_sq;

  // Calculate also the change in the momentum squared also ???
  // G4double veloc_square = y_posvel.GetVelocity().mag2();
  // ...

  if( dyerr_pos_sq > ( dyerr_mom_rel_sq * hstep*hstep ) )
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
G4double GXMagInt_Driver_ComputeNewStepSize(GXMagInt_Driver *This, 
                          G4double  errMaxNorm,    // max error  (normalised)
                          G4double  hstepCurrent)  // current step size
{
  G4double hnew;

  // Compute size of next Step for a failed step
  if(errMaxNorm > 1.0 )
  {
    // Step failed; compute the size of retrial Step.
    hnew = GXMagInt_Driver_GetSafety(This)*hstepCurrent*
      pow(errMaxNorm,GXMagInt_Driver_GetPshrnk(This)) ;
  } else if(errMaxNorm > 0.0 ) {
    // Compute size of next Step for a successful step
    hnew = GXMagInt_Driver_GetSafety(This)*hstepCurrent*
      pow(errMaxNorm,GXMagInt_Driver_GetPgrow(This)) ;
  } else {
    // if error estimate is zero (possible) or negative (dubious)
    hnew = max_stepping_increase * hstepCurrent; 
  }

  return hnew;
}

// ---------------------------------------------------------------------------
//G4MagIntegratorDriver.icc

FQUALIFIER
G4double GXMagInt_Driver_GetHmin( GXMagInt_Driver *This )
{
  return This->fMinimumStep;
} 

FQUALIFIER
G4double GXMagInt_Driver_GetSafety( GXMagInt_Driver *This )
{
  return This->safety;
}

FQUALIFIER
G4double GXMagInt_Driver_GetPshrnk( GXMagInt_Driver *This )
{
  return This->pshrnk;
} 

FQUALIFIER
G4double GXMagInt_Driver_GetPgrow( GXMagInt_Driver *This )
{
  return This->pgrow;
}
 
FQUALIFIER
G4double GXMagInt_Driver_GetErrcon( GXMagInt_Driver *This )
{
  return This->errcon;
}

FQUALIFIER
const GXClassicalRK4* GXMagInt_Driver_GetStepper(GXMagInt_Driver *This)
{
  return This->pIntStepper;
}

FQUALIFIER
void GXMagInt_Driver_GetDerivatives( GXMagInt_Driver *This,
				     GXFieldTrack &y_curr,// const, INput
				     G4double     dydx[])  // OUTput
{ 
  G4double  tmpValArr[GXFieldTrack_ncompSVEC];
  GXFieldTrack_DumpToArray(&y_curr, tmpValArr );
  GXClassicalRK4_RightHandSide(This->pIntStepper, tmpValArr , dydx );
}

//@@@G4FWP test implemenation and should be used only for performance evaluation
FQUALIFIER
G4bool GXMagInt_Driver_AccurateAdvance2( GXMagInt_Driver *This,
					 GXFieldTrack& y_current,
					 G4double     hstep,
					 G4double     eps,
					 G4double hinitial )
{
  // Runge-Kutta driver with adaptive stepsize control. Integrate starting
  // values at y_current over hstep x2 with accuracy eps. 
  // On output ystart is replaced by values at the end of the integration 
  // interval. RightHandSide is the right-hand side of ODE system. 
  // The source is similar to odeint routine from NRC p.721-722 .

  //  G4int nstp;
  G4int i =0;
  G4double x, hnext, hdid, h;

  G4double y[GXFieldTrack_ncompSVEC], dydx[GXFieldTrack_ncompSVEC];
  //  G4double ystart[GXFieldTrack_ncompSVEC], yEnd[GXFieldTrack_ncompSVEC]; 
  G4double ystart[GXFieldTrack_ncompSVEC];
  G4double x1;
  //  G4double x2;
  G4bool succeeded = true;

  G4double startCurveLength;

  //  Ensure that hstep > 0
  //
  if( hstep <= 0.0 )
  { 
    if(hstep==0.0)
    {
      //warning message
      return succeeded; 
    }
    else
    { 
      //exception
      return false;
    }
  }

  GXFieldTrack_DumpToArray(&y_current, ystart );

  //  startCurveLength= y_current.GetCurveLength();
  startCurveLength= GXFieldTrack_GetCurveLength(&y_current);
  x1= startCurveLength; 
  //  x2= x1 + hstep;

  if ( (hinitial > 0.0) && (hinitial < hstep)
    && (hinitial > perMillion * hstep) )
  {
     h = hinitial;
  }
  else  //  Initial Step size "h" defaults to the full interval
  {
     h = hstep;
  }

  //  h = hstep;

  x = x1;

  for (i=0 ;i < GXFieldTrack_ncompSVEC ;i++)  { y[i] = ystart[i]; }

  //  G4bool lastStep= false;
  //  nstp=1;

  //  do
  {
    GPThreeVector StartPos = GPThreeVector_create( y[0], y[1], y[2] );

    GXClassicalRK4_RightHandSide(This->pIntStepper, y, dydx );

    // Perform the Integration
    //      
    if( h > This->fMinimumStep )
    { 
      GXMagInt_Driver_OneGoodStep2(This, y,dydx,x,h,eps,hdid,hnext) ;
    }
    else
    {
      GXFieldTrack yFldTrk;
      GXFieldTrack_Constructor(&yFldTrk, 
			       GPThreeVector_create(0,0,0), 
			       GPThreeVector_create(0,0,0), 
			       0., 0., 0., 0.);

      //      G4double dchord_step, dyerr, dyerr_len;   // What to do with these ?
      G4double dchord_step, dyerr_len;   // What to do with these ?
      GXFieldTrack_LoadFromArray(&yFldTrk, y); 
      GXFieldTrack_SetCurveLength(&yFldTrk, x );

      GXMagInt_Driver_QuickAdvance(This, 
				   yFldTrk, dydx, h, dchord_step, dyerr_len ); 
      //-----------------------------------------------------

      GXFieldTrack_DumpToArray(&yFldTrk,y);    

      /*
      if( h == 0.0 )
      { 
	; //invoke exception
      }

      dyerr = dyerr_len / h;
      hdid= h;
      x += hdid;

      // Compute suggested new step
      hnext= GXMagInt_Driver_ComputeNewStepSize(This, dyerr/eps, h);
      */
    }

    //  Avoid numerous small last steps
    /*
    if( (h < eps * hstep) || (h < This->fSmallestFraction * startCurveLength) )
    {
      // No more integration -- the next step will not happen
      lastStep = true;  
    }
    else
    {
      // Check the proposed next stepsize
      if(fabs(hnext) <= GXMagInt_Driver_GetHmin(This))
      {
        // Make sure that the next step is at least Hmin.
        h = GXMagInt_Driver_GetHmin(This);
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
    */
  } 
  //while ( ((nstp++)<=This->fMaxNoSteps) && (x < x2) && (!lastStep) );
     // Have we reached the end ?
     // --> a better test might be x-x2 > an_epsilon

  /*
  succeeded=  (x>=x2);  // If it was a "forced" last step

  for (i=0; i < GXFieldTrack_ncompSVEC ;i++)  { yEnd[i] = y[i]; }

  // Put back the values.
  GXFieldTrack_LoadFromArray(&y_current, yEnd );
  GXFieldTrack_SetCurveLength(&y_current, x );

  if(nstp > This->fMaxNoSteps)
  {
    succeeded = false;
  }
  */

  return succeeded;
}  // end of AccurateAdvance ...........................
