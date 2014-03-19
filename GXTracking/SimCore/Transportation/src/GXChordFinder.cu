#include "GXChordFinder.h"
#include "GXMagneticField.h"
#include "GPConstants.h"
#include "GPUtils.h"
#include "stdio.h"

FQUALIFIER
void GXChordFinder_Constructor(GXChordFinder *This,
			       GXMagInt_Driver* pIntegrationDriver)
{
  This->fDeltaChord = 0.25 * millimeter ; // default delta Chord 
  This->fFirstFraction = 0.999; 
  This->fFractionLast = 1.00;
  This->fFractionNextEstimate = 0.98;
  This->fLastStepEstimate_Unconstrained = DBL_MAX;

  This->fEquation = 0;
  This->fDriversStepper = 0;
  This->fIntgrDriver= pIntegrationDriver;

}

// ..........................................................................
FQUALIFIER
void GXChordFinder_Constructor2( GXChordFinder     *This,
				 GXMagneticField*  theMagField,
				 G4double          stepMinimum, 
				 GXClassicalRK4*   pItsStepper )
{
  This->fDeltaChord = 0.25 * millimeter ; // default delta Chord
  This->fFirstFraction = 0.999; 
  This->fFractionLast = 1.00;  
  This->fFractionNextEstimate = 0.98; 
  This->fLastStepEstimate_Unconstrained = DBL_MAX; 

  This->fDriversStepper =0;                    // Dependent objects 
  This->fEquation = 0;      

  //  Construct the Chord Finder
  //  by creating in inverse order the  Driver, the Stepper and EqRhs ...

  GXEquationOfMotion pEquation;
  GXEquationOfMotion_Constructor(&pEquation,theMagField,-1.0); 
  This->fEquation = &pEquation;                            

  if( pItsStepper == 0 )
  { 
    GXClassicalRK4_Constructor(This->fDriversStepper, &pEquation);
    pItsStepper = This->fDriversStepper;
  }

  GXMagInt_Driver_Constructor(This->fIntgrDriver, stepMinimum, pItsStepper);

}

// this is the master method that is called from GXPropagatorInField

FQUALIFIER
G4double 
GXChordFinder_AdvanceChordLimited( GXChordFinder *This,
				   GXFieldTrack& yCurrent,
				   G4double      stepMax,
				   G4double      epsStep) 
                                   //rest arguments are not actually used
{
  G4double stepPossible;
  G4double dyErr;
  GXFieldTrack&  yEnd =  yCurrent;
  G4double  startCurveLen= GXFieldTrack_GetCurveLength(&yCurrent);
  G4double nextStep;

  //*****************************************************************
  stepPossible= GXChordFinder_FindNextChord(This, yCurrent, stepMax, 
					    yEnd, dyErr, epsStep,
					    &nextStep);
  //*****************************************************************

  G4bool good_advance;

  if ( dyErr < epsStep * stepPossible )
  {
     // Accept this accuracy.
     yCurrent = yEnd;
     good_advance = true; 
  }
  else
  {  
    // Advance more accurately to "end of chord"
    //*****************************************************************
    good_advance = GXMagInt_Driver_AccurateAdvance(This->fIntgrDriver, 
						   yCurrent, stepPossible,
						   epsStep, nextStep);
    //*****************************************************************
    if ( ! good_advance )
    { 
      // In this case the driver could not do the full distance
      stepPossible= GXFieldTrack_GetCurveLength(&yCurrent)-startCurveLen;
    }
  }
  return stepPossible;
}

//@@@G4FWP test implemenation and should be used only for performance evaluation
FQUALIFIER
G4double 
GXChordFinder_AdvanceChordLimited2( GXChordFinder *This,
				    GXFieldTrack& yCurrent,
				    G4double      stepMax,
				    G4double      epsStep) 
{ 

  G4double stepPossible;

  G4double dyErr;
  GXFieldTrack&  yEnd =  yCurrent;
  G4double  startCurveLen= GXFieldTrack_GetCurveLength(&yCurrent);

  G4double nextStep;

  //*****************************************************************
  stepPossible= GXChordFinder_FindNextChord2(This, yCurrent, stepMax, 
					     yEnd, dyErr, epsStep,
					     &nextStep);
  //*****************************************************************

  G4bool good_advance;

  nextStep = 0.0;

  if ( dyErr < epsStep * stepPossible )
  {
     // Accept this accuracy.
     yCurrent = yEnd;
     good_advance = true; 
  }
  else
  {  
    // Advance more accurately to "end of chord"
    //*****************************************************************
    good_advance = GXMagInt_Driver_AccurateAdvance2(This->fIntgrDriver, 
						    yCurrent, stepPossible,
						    epsStep, nextStep);
    //*****************************************************************
    if ( ! good_advance )
      { 
      // In this case the driver could not do the full distance
	stepPossible= GXFieldTrack_GetCurveLength(&yCurrent)-startCurveLen;
      }
  }
  return stepPossible;
}

// ............................................................................

FQUALIFIER
G4double
GXChordFinder_FindNextChord( GXChordFinder *This, 
			     const  GXFieldTrack& yStart,
			     G4double     stepMax,
			     GXFieldTrack&   yEnd, // Endpoint
			     G4double&   dyErrPos, // Error of endpoint
			     G4double    epsStep,
			     G4double*  pStepForAccuracy)
{
  // Returns Length of Step taken

  GXFieldTrack yCurrent=  yStart;  
  G4double    stepTrial, stepForAccuracy;
  G4double    dydx[GXFieldTrack_ncompSVEC]; 

  //  1.)  Try to "leap" to end of interval
  //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
  // 2a.)  If d_chord is not good enough, find one that is.
  
  G4bool    validEndPoint= false;
  G4double  dChordStep, lastStepLength; //  stepOfLastGoodChord;

  GXMagInt_Driver_GetDerivatives(This->fIntgrDriver, yCurrent, dydx );

  G4int     noTrials=0;
  const G4double safetyFactor= This->fFirstFraction; //0.975 or 0.99 ? was 0.999

  stepTrial = GPfmin( stepMax, 
		     safetyFactor*This->fLastStepEstimate_Unconstrained );

  G4double newStepEst_Uncons= 0.0; 
  do
  { 
     G4double stepForChord;  
     yCurrent = yStart;    // Always start from initial point
    
     //            ************
     GXMagInt_Driver_QuickAdvance(This->fIntgrDriver, yCurrent, dydx, 
				  stepTrial, dChordStep, dyErrPos);
     //            ************
     
     // We check whether the criterion is met here.
     // dChordStep <= fDeltaChord = fDefaultDeltaChord = 0.25*millimeter  
     validEndPoint = GXChordFinder_AcceptableMissDist(This,dChordStep);

     lastStepLength = stepTrial; 

     // This method estimates to step size for a good chord.
     stepForChord = GXChordFinder_NewStep(This, stepTrial, 
					  dChordStep, newStepEst_Uncons );

     if( ! validEndPoint )
     {
        if( stepTrial<=0.0 )
        {
          stepTrial = stepForChord;
        }
        else if (stepForChord <= stepTrial)
        {
          // Reduce by a fraction, possibly up to 20% 
          stepTrial = GPfmin( stepForChord, This->fFractionLast * stepTrial);
        }
        else
        {
          stepTrial *= 0.1;
        }
     }
     if (noTrials > 100) {
        printf("GPU-ERROR: too many trials (100) in GXChordFinder_FindNextChord %e %e\n",
               stepTrial,stepMax);
        break;
     }
     noTrials++; 
  }
  while( ! validEndPoint );   // End of do-while  RKD 

  if( newStepEst_Uncons > 0.0  )
  {
     This->fLastStepEstimate_Unconstrained= newStepEst_Uncons;
  }

  if( pStepForAccuracy )
  { 
     // Calculate the step size required for accuracy, if it is needed
     //
     G4double dyErr_relative = dyErrPos/(epsStep*lastStepLength);
     if( dyErr_relative > 1.0 )
     {
       stepForAccuracy = 
	 GXMagInt_Driver_ComputeNewStepSize(This->fIntgrDriver,
					    dyErr_relative,
					    lastStepLength );
     }
     else
     {
        stepForAccuracy = 0.0;   // Convention to show step was ok 
     }
     *pStepForAccuracy = stepForAccuracy;
  }

  yEnd=  yCurrent;  
  return stepTrial; 
}

//@@@G4FWP test implemenation and should be used only for performance evaluation
FQUALIFIER
G4double
GXChordFinder_FindNextChord2( GXChordFinder *This, 
			      const  GXFieldTrack& yStart,
			      G4double     stepMax,
			      GXFieldTrack&   yEnd, // Endpoint
			      G4double&   dyErrPos, // Error of endpoint
			      G4double    epsStep,
			      G4double*  pStepForAccuracy)
{
  // Returns Length of Step taken

  GXFieldTrack yCurrent=  yStart;  
  G4double    stepTrial, stepForAccuracy;
  G4double    dydx[GXFieldTrack_ncompSVEC]; 

  //  1.)  Try to "leap" to end of interval
  //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
  // 2a.)  If d_chord is not good enough, find one that is.
  
  //  G4bool    validEndPoint= false;
  G4double  dChordStep, lastStepLength; //  stepOfLastGoodChord;

  GXMagInt_Driver_GetDerivatives(This->fIntgrDriver, yCurrent, dydx );

  G4int     noTrials=0;
  const G4double safetyFactor= This->fFirstFraction; //0.975 or 0.99 ? was 0.999

  stepTrial = GPfmin( stepMax, 
		     safetyFactor*This->fLastStepEstimate_Unconstrained );

  G4double newStepEst_Uncons= 0.0; 
  //  do
  { 
    //     G4double stepForChord;  
     yCurrent = yStart;    // Always start from initial point
    
     //            ************
     GXMagInt_Driver_QuickAdvance(This->fIntgrDriver, yCurrent, dydx, 
				  stepTrial, dChordStep, dyErrPos);
     //            ************
     
     // We check whether the criterion is met here.
     // dChordStep <= fDeltaChord = fDefaultDeltaChord = 0.25*millimeter  
     //     validEndPoint = GXChordFinder_AcceptableMissDist(This,dChordStep);

     lastStepLength = stepTrial; 

     // This method estimates to step size for a good chord.
     //     stepForChord = GXChordFinder_NewStep(This, stepTrial, 
     //					  dChordStep, newStepEst_Uncons );

     /*
     if( ! validEndPoint )
     {
        if( stepTrial<=0.0 )
        {
          stepTrial = stepForChord;
        }
        else if (stepForChord <= stepTrial)
        {
          // Reduce by a fraction, possibly up to 20% 
          stepTrial = GPfmin( stepForChord, This->fFractionLast * stepTrial);
        }
        else
        {
          stepTrial *= 0.1;
        }
     }
     */
     noTrials++; 
  }
  //  while( ! validEndPoint );   // End of do-while  RKD 

  if( newStepEst_Uncons > 0.0  )
  {
     This->fLastStepEstimate_Unconstrained= newStepEst_Uncons;
  }

 if( pStepForAccuracy )
  { 
     // Calculate the step size required for accuracy, if it is needed
     //
     G4double dyErr_relative = dyErrPos/(epsStep*lastStepLength);
     if( dyErr_relative > 1.0 )
     {
       stepForAccuracy = 
	 GXMagInt_Driver_ComputeNewStepSize(This->fIntgrDriver,
					    dyErr_relative,
					    lastStepLength );
     }
     else
     {
        stepForAccuracy = 0.0;   // Convention to show step was ok 
     }
     *pStepForAccuracy = stepForAccuracy;
  }

  yEnd=  yCurrent;  
  return stepTrial; 
}

// ...........................................................................

FQUALIFIER
G4double GXChordFinder_NewStep( GXChordFinder *This,
				G4double  stepTrialOld, 
                                G4double  dChordStep, // Curr. dchord achieved
                                G4double& stepEstimate_Unconstrained )  
{
  // Is called to estimate the next step size, even for successful steps,
  // in order to predict an accurate 'chord-sensitive' first step
  // which is likely to assist in more performant 'stepping'.

  G4double stepTrial;

  if (dChordStep > 0.0)
  {
    stepEstimate_Unconstrained =
                 stepTrialOld*sqrt( This->fDeltaChord / dChordStep );
    stepTrial =  This->fFractionNextEstimate * stepEstimate_Unconstrained;
  }
  else
  {
    // Should not update the Unconstrained Step estimate: incorrect!
    stepTrial =  stepTrialOld * 2.; 
  }

  if( stepTrial <= 0.001 * stepTrialOld)
  {
     if ( dChordStep > 1000.0 * This->fDeltaChord )
     {
        stepTrial= stepTrialOld * 0.03;   
     }
     else
     {
        if ( dChordStep > 100. * This->fDeltaChord )
        {
          stepTrial= stepTrialOld * 0.1;   
        }
        else   // Try halving the length until dChordStep OK
        {
          stepTrial= stepTrialOld * 0.5;   
        }
     }
  }
  else if (stepTrial > 1000.0 * stepTrialOld)
  {
     stepTrial= 1000.0 * stepTrialOld;
  }

  if( stepTrial == 0.0 )
  {
     stepTrial= 0.000001;
  }
  return stepTrial;
}


// ...........................................................................

FQUALIFIER
GXFieldTrack
GXChordFinder_ApproxCurvePointS( GXChordFinder *This, 
				 const GXFieldTrack&  CurveA_PointVelocity, 
				 const GXFieldTrack&  CurveB_PointVelocity, 
				 const GXFieldTrack&  ApproxCurveV,
				 const GPThreeVector& CurrentE_Point,
				 const GPThreeVector& CurrentF_Point,
				 const GPThreeVector& PointG,
				 G4bool first, G4double eps_step)
{
  // ApproxCurvePointS is 2nd implementation of ApproxCurvePoint.
  // Use Brent Algorithm (or InvParabolic) when possible.
  // Given a starting curve point A (CurveA_PointVelocity), curve point B
  // (CurveB_PointVelocity), a point E which is (generally) not on the curve
  // and  a point F which is on the curve (first approximation), find new
  // point S on the curve closer to point E. 
  // While advancing towards S utilise 'eps_step' as a measure of the
  // relative accuracy of each Step.

  //  GXFieldTrack EndPoint(CurveA_PointVelocity);
  GXFieldTrack EndPoint = CurveA_PointVelocity;
  if(!first){EndPoint= ApproxCurveV;}

  GPThreeVector Point_A,Point_B;
  Point_A=GXFieldTrack_GetPosition(&CurveA_PointVelocity);
  Point_B=GXFieldTrack_GetPosition(&CurveB_PointVelocity);

  G4double xa,xb,xc,ya,yb,yc;
 
  // InverseParabolic. AF Intersects (First Part of Curve) 

  if(first)
  {
    xa=0.;
    ya=GPThreeVector_mag(GPThreeVector_sub(PointG,Point_A));
    xb=GPThreeVector_mag(GPThreeVector_sub(Point_A,CurrentF_Point));
    yb=-GPThreeVector_mag(GPThreeVector_sub(PointG,CurrentF_Point));
    xc=GPThreeVector_mag(GPThreeVector_sub(Point_A,Point_B));
    yc=-1.0*GPThreeVector_mag(GPThreeVector_sub(CurrentE_Point,Point_B));
  }    
  else
  {
    xa=0.;
    ya=GPThreeVector_mag(GPThreeVector_sub(Point_A,CurrentE_Point));
    xb=GPThreeVector_mag(GPThreeVector_sub(Point_A,CurrentF_Point));
    yb=GPThreeVector_mag(GPThreeVector_sub(PointG,CurrentF_Point));
    xc=GPThreeVector_mag(GPThreeVector_sub(Point_A,Point_B));
    yc=-1.0*GPThreeVector_mag(GPThreeVector_sub(Point_B,PointG));
    if(xb==0.)
    {
      EndPoint=
	GXChordFinder_ApproxCurvePointV(This, 
					CurveA_PointVelocity, 
					CurveB_PointVelocity,
					CurrentE_Point, eps_step);
      return EndPoint;
    }
  }

  const G4double tolerance= 1.e-12;
  if( abs(ya)<=tolerance|| abs(yc)<=tolerance)
  {
    ; // What to do for the moment: return the same point as at start
      // then PropagatorInField will take care
  }
  else
  {
    G4double test_step = GXChordFinder_InvParabolic(xa,ya,xb,yb,xc,yc);
    G4double curve;
    if(first)
    {
      curve=abs(GXFieldTrack_GetCurveLength(&EndPoint)
                    -GXFieldTrack_GetCurveLength(&ApproxCurveV));
    }
    else
    {
      test_step=(test_step-xb);
      curve=abs(GXFieldTrack_GetCurveLength(&EndPoint)
                    -GXFieldTrack_GetCurveLength(&CurveB_PointVelocity));
      xb=GPThreeVector_mag(GPThreeVector_sub(CurrentF_Point,Point_B));
    }
      
    if(test_step<=0)    { test_step=0.1*xb; }
    if(test_step>=xb)   { test_step=0.5*xb; }
    if(test_step>=curve){ test_step=0.5*curve; } 

    if(curve*(1.+eps_step)<xb) // Similar to ReEstimate Step from
    {                          // G4VIntersectionLocator
      test_step=0.5*curve;
    }

    GXMagInt_Driver_AccurateAdvance(This->fIntgrDriver,
				    EndPoint,test_step, eps_step, 0);
      
  }
  return EndPoint;
}


// ...........................................................................

FQUALIFIER
GXFieldTrack 
GXChordFinder_ApproxCurvePointV( GXChordFinder *This, 
				 const GXFieldTrack& CurveA_PointVelocity, 
				 const GXFieldTrack& CurveB_PointVelocity, 
				 const GPThreeVector& CurrentE_Point,
				 G4double eps_step)
{
  // If r=|AE|/|AB|, and s=true path lenght (AB)
  // return the point that is r*s along the curve!
 
  GXFieldTrack   Current_PointVelocity = CurveA_PointVelocity; 

  GPThreeVector  CurveA_Point= GXFieldTrack_GetPosition(&CurveA_PointVelocity);
  GPThreeVector  CurveB_Point= GXFieldTrack_GetPosition(&CurveB_PointVelocity);

  GPThreeVector  ChordAB_Vector= GPThreeVector_sub(CurveB_Point,CurveA_Point);
  GPThreeVector  ChordAE_Vector= GPThreeVector_sub(CurrentE_Point,CurveA_Point);

  G4double       ABdist= GPThreeVector_mag(ChordAB_Vector);
  G4double  curve_length;  //  A curve length  of AB
  G4double  AE_fraction; 
  
  curve_length= GXFieldTrack_GetCurveLength(&CurveB_PointVelocity)
              - GXFieldTrack_GetCurveLength(&CurveA_PointVelocity);  
 
  G4double  integrationInaccuracyLimit= GPfmax( perMillion, 0.5*eps_step ); 
  if( curve_length < ABdist * (1. - integrationInaccuracyLimit) )
  { 
    // Take default corrective action: adjust the maximum curve length. 
    // NOTE: this case only happens for relatively straight paths.
    // curve_length = ABdist; 
  }

  G4double  new_st_length; 

  if ( ABdist > 0.0 )
  {
    AE_fraction = GPThreeVector_mag(ChordAE_Vector) / ABdist;
  }
  else
  {
     AE_fraction = 0.5;                         // Guess .. ?; 
  }
  
  if( (AE_fraction> 1.0 + perMillion) || (AE_fraction< 0.) )
  {
     // This course can now result if B has been re-evaluated, 
     // without E being recomputed (1 July 99).
     // In this case this is not a "real error" - but it is undesired
     // and we cope with it by a default corrective action ...
     //
     AE_fraction = 0.5;                         // Default value
  }

  new_st_length= AE_fraction * curve_length; 

  if ( AE_fraction > 0.0 )
  { 
    GXMagInt_Driver_AccurateAdvance(This->fIntgrDriver,
				    Current_PointVelocity, 
				    new_st_length, eps_step ,0);
     //
     // In this case it does not matter if it cannot advance the full distance
  }

  // If there was a memory of the step_length actually required at the start 
  // of the integration Step, this could be re-used ...

  //  G4cout.precision(14);

  return Current_PointVelocity;
}

FQUALIFIER
GXMagInt_Driver* GXChordFinder_GetIntegrationDriver(GXChordFinder *This)
{
  return This->fIntgrDriver;
}

FQUALIFIER
G4bool GXChordFinder_AcceptableMissDist(GXChordFinder *This,
					G4double dChordStep)
{ 
  return (dChordStep <= This->fDeltaChord) ;
}

FQUALIFIER
G4double GXChordFinder_GetDeltaChord(GXChordFinder *This)
{
  return This->fDeltaChord;
}

FQUALIFIER
void GXChordFinder_SetDeltaChord(GXChordFinder *This,G4double newval)
{
  This->fDeltaChord=newval;
}
// ......................................................................

FQUALIFIER 
G4double GXChordFinder_GetLastStepEstimateUnc(GXChordFinder *This)
{
  return This->fLastStepEstimate_Unconstrained;   
} 

FQUALIFIER 
void GXChordFinder_SetLastStepEstimateUnc(GXChordFinder *This,
					  G4double stepEst )
{
  This->fLastStepEstimate_Unconstrained = stepEst;    
}

FQUALIFIER
void GXChordFinder_ResetStepEstimate(GXChordFinder *This)
{
  This->fLastStepEstimate_Unconstrained = DBL_MAX;    
}

// ......................................................................

FQUALIFIER 
G4double GXChordFinder_GetFirstFraction(GXChordFinder *This) 
{ return This->fFirstFraction; } 

FQUALIFIER 
G4double GXChordFinder_GetFractionLast(GXChordFinder *This)  
{ return This->fFractionLast; } 

FQUALIFIER 
void GXChordFinder_SetFirstFraction(GXChordFinder *This,
				    G4double val)
{ This->fFirstFraction=val; }

FQUALIFIER 
G4double GXChordFinder_InvParabolic ( const G4double xa, 
				      const G4double ya,
				      const G4double xb, 
				      const G4double yb,
				      const G4double xc, 
				      const G4double yc )
{       
  const G4double R = yb/yc;
  const G4double S = yb/ya;
  const G4double T = ya/yc;
  const G4double Q = (T-1)*(R-1)*(S-1);
  if (fabs(Q) < DBL_MIN ) return  DBL_MAX;
  
  const G4double P = S*(T*(R-T)*(xc-xb) - (1-R)*(xb-xa));
  return xb + P/Q;
}
