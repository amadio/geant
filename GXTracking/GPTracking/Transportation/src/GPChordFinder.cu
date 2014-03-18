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
// $Id: G4ChordFinder.cc,v 1.53 2009-05-18 14:22:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// 25.02.97 - John Apostolakis - Design and implementation 
// -------------------------------------------------------------------

//#include <iomanip>

#include "GPChordFinder.h"
#include "GPMagneticField.h"
#include "GPConstants.h"
#include "GPUtils.h"

//#include "G4Mag_UsualEqRhs.hh"
//#include "GPEquationOfMotion.h"
//#include "GPClassicalRK4.h"
// ..........................................................................

FQUALIFIER
void GPChordFinder_Constructor(GPChordFinder *This,
			       GPMagInt_Driver* pIntegrationDriver)
{
  This->fDefaultDeltaChord = 0.25 * millimeter ;    // Parameters
  This->fDeltaChord = This->fDefaultDeltaChord ;    //   Internal parameters
  This->fFirstFraction = 0.999; 
  This->fFractionLast = 1.00;  
  This->fFractionNextEstimate = 0.98; 
  This->fMultipleRadius = 15.0;
  This->fStatsVerbose =0;
  This->fDriversStepper =0;                    // Dependent objects 
  This->fAllocatedStepper = false;
  This->fEquation = 0;      
  This->fTotalNoTrials_FNC = 0;
  This->fNoCalls_FNC =0; 
  This->fmaxTrials_FNC =0;

  // Simple constructor -- it does not create equation
  This->fIntgrDriver= pIntegrationDriver;
  This->fAllocatedStepper= false;

  This->fLastStepEstimate_Unconstrained = DBL_MAX;  // Should move q, p to

  GPChordFinder_SetFractions_Last_Next(This,
				       This->fFractionLast, 
				       This->fFractionNextEstimate);  
  // check the values and set the other parameters
}


// ..........................................................................
//			    G4MagIntegratorStepper* pItsStepper )

FQUALIFIER
void GPChordFinder_Constructor2( GPChordFinder *This,
				 GPMagneticField*        theMagField,
				 G4double                stepMinimum, 
				 GPClassicalRK4* pItsStepper )
{
  This->fDefaultDeltaChord = 0.25 * millimeter ;    // Constants
  This->fDeltaChord = This->fDefaultDeltaChord ;    // Parameters
  This->fFirstFraction = 0.999; 
  This->fFractionLast = 1.00;  
  This->fFractionNextEstimate = 0.98; 
  This->fMultipleRadius = 15.0;
  This->fStatsVerbose =0;
  This->fDriversStepper =0;                    // Dependent objects 
  This->fAllocatedStepper = false;
  This->fEquation = 0;      
  This->fTotalNoTrials_FNC = 0;  // State - stats
  This->fNoCalls_FNC =0; 
  This->fmaxTrials_FNC =0; 

  //  Construct the Chord Finder
  //  by creating in inverse order the  Driver, the Stepper and EqRhs ...

  //  G4Mag_EqRhs *pEquation = new G4Mag_UsualEqRhs(theMagField);
  GPEquationOfMotion pEquation;
  GPEquationOfMotion_Constructor(&pEquation,theMagField); 

  //  fEquation = pEquation;                            
  This->fEquation = &pEquation;                            
  This->fLastStepEstimate_Unconstrained = DBL_MAX;   // Should move q, p to
                                                     //    G4FieldTrack ??

  GPChordFinder_SetFractions_Last_Next( This,
					This->fFractionLast, 
					This->fFractionNextEstimate);  
    // check the values and set the other parameters

  // --->>  Charge    Q = 0 
  // --->>  Momentum  P = 1       NOMINAL VALUES !!!!!!!!!!!!!!!!!!

  if( pItsStepper == 0 )
  { 
    //     pItsStepper = fDriversStepper = new G4ClassicalRK4(pEquation);
    GPClassicalRK4_Constructor(This->fDriversStepper, &pEquation, 0);
    pItsStepper = This->fDriversStepper;
    This->fAllocatedStepper= true;
  }
  else
  {
     This->fAllocatedStepper= false; 
  }
  //  fIntgrDriver = new GPMagInt_Driver(stepMinimum, pItsStepper, 
  //                                     pItsStepper->GetNumberOfVariables() );
  GPMagInt_Driver_Constructor(This->fIntgrDriver, stepMinimum, pItsStepper, 
			      GPClassicalRK4_GetNumberOfVariables(pItsStepper),
			      0);
}


// ......................................................................

//GPChordFinder::~GPChordFinder()
//{
//  delete   fEquation; // fIntgrDriver->pIntStepper->theEquation_Rhs;
//  if( fAllocatedStepper)
//  { 
//     delete fDriversStepper; 
//  }
//  delete   fIntgrDriver; 
//
//  if( fStatsVerbose ) { PrintStatistics(); }
//}


// ......................................................................

FQUALIFIER
void   
GPChordFinder_SetFractions_Last_Next( GPChordFinder *This, 
				      G4double fractLast, G4double fractNext )
{ 
  // Use -1.0 as request for Default.
  if( fractLast == -1.0 )   fractLast = 1.0;   // 0.9;
  if( fractNext == -1.0 )   fractNext = 0.98;  // 0.9; 

  // fFirstFraction  = 0.999; // Orig 0.999 A safe value, range: ~ 0.95 - 0.999
  // fMultipleRadius = 15.0;  // For later use, range: ~  2 - 20 

  if( This->fStatsVerbose )
  { 
    //    G4cout << " ChordFnd> Trying to set fractions: "
    //           << " first " << fFirstFraction
    //           << " last " <<  fractLast
    //           << " next " <<  fractNext
    //           << " and multiple " << fMultipleRadius
    //           << G4endl;
  } 

  if( (fractLast > 0.0) && (fractLast <=1.0) ) 
  {
    This->fFractionLast= fractLast;
  }
  else
  {
    //    G4cerr << "GPChordFinder::SetFractions_Last_Next: Invalid "
    //           << " fraction Last = " << fractLast
    //           << " must be  0 <  fractionLast <= 1 " << G4endl;
  }
  if( (fractNext > 0.0) && (fractNext <1.0) )
  {
    This->fFractionNextEstimate = fractNext;
  }
  else
  {
    ;
    //    G4cerr << "GPChordFinder:: SetFractions_Last_Next: Invalid "
    //           << " fraction Next = " << fractNext
    //           << " must be  0 <  fractionNext < 1 " << G4endl;
  }
}


// ......................................................................

FQUALIFIER
G4double 
GPChordFinder_AdvanceChordLimited( GPChordFinder *This,
				   GPFieldTrack& yCurrent,
				   G4double      stepMax,
				   G4double      epsStep,
				   const GPThreeVector latestSafetyOrigin,
				   G4double       latestSafetyRadius )
{
  G4double stepPossible;
  G4double dyErr;
  //  G4FieldTrack yEnd( yCurrent);
  //  GPFieldTrack& yEnd;
  GPFieldTrack&  yEnd =  yCurrent;
  G4double  startCurveLen= GPFieldTrack_GetCurveLength(&yCurrent);
  G4double nextStep;
  //            *************
  stepPossible= GPChordFinder_FindNextChord(This, yCurrent, stepMax, yEnd, 
					    dyErr, epsStep,
                              &nextStep, latestSafetyOrigin, latestSafetyRadius
                             );
  //            *************

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
     //                           ***************
    good_advance = GPMagInt_Driver_AccurateAdvance(This->fIntgrDriver, 
						   yCurrent, stepPossible,
						   epsStep, nextStep);
    if ( ! good_advance )
    { 
      // In this case the driver could not do the full distance
      stepPossible= GPFieldTrack_GetCurveLength(&yCurrent)-startCurveLen;
    }
  }
  return stepPossible;
}


// ............................................................................

FQUALIFIER
G4double
GPChordFinder_FindNextChord( GPChordFinder *This, 
			     const  GPFieldTrack& yStart,
			     G4double     stepMax,
			     GPFieldTrack&   yEnd, // Endpoint
			     G4double&   dyErrPos, // Error of endpoint
			     G4double    epsStep,
			     G4double*  pStepForAccuracy, 
			     const  GPThreeVector, //  latestSafetyOrigin,
			     G4double       //  latestSafetyRadius 
			     )
{
  // Returns Length of Step taken

  GPFieldTrack yCurrent=  yStart;  
  G4double    stepTrial, stepForAccuracy;
  G4double    dydx[GPFieldTrack_ncompSVEC]; 

  //  1.)  Try to "leap" to end of interval
  //  2.)  Evaluate if resulting chord gives d_chord that is good enough.
  // 2a.)  If d_chord is not good enough, find one that is.
  
  G4bool    validEndPoint= false;
  G4double  dChordStep, lastStepLength; //  stepOfLastGoodChord;

  GPMagInt_Driver_GetDerivatives(This->fIntgrDriver, yCurrent, dydx );

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
     GPMagInt_Driver_QuickAdvance(This->fIntgrDriver, yCurrent, dydx, 
				  stepTrial, dChordStep, dyErrPos);
     //            ************
     
     // We check whether the criterion is met here.
     // dChordStep <= fDeltaChord = fDefaultDeltaChord = 0.25*millimeter  
     validEndPoint = GPChordFinder_AcceptableMissDist(This,dChordStep);

     lastStepLength = stepTrial; 

     // This method estimates to step size for a good chord.
     stepForChord = GPChordFinder_NewStep(This, stepTrial, 
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
     noTrials++; 
  }
  while( ! validEndPoint );   // End of do-while  RKD 

  if( newStepEst_Uncons > 0.0  )
  {
     This->fLastStepEstimate_Unconstrained= newStepEst_Uncons;
  }

  GPChordFinder_AccumulateStatistics(This, noTrials );

  if( pStepForAccuracy )
  { 
     // Calculate the step size required for accuracy, if it is needed
     //
     G4double dyErr_relative = dyErrPos/(epsStep*lastStepLength);
     if( dyErr_relative > 1.0 )
     {
       stepForAccuracy = 
	 GPMagInt_Driver_ComputeNewStepSize(This->fIntgrDriver,
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
G4double GPChordFinder_NewStep( GPChordFinder *This,
				G4double  stepTrialOld, 
                                G4double  dChordStep, // Curr. dchord achieved
                                G4double& stepEstimate_Unconstrained )  
{
  // Is called to estimate the next step size, even for successful steps,
  // in order to predict an accurate 'chord-sensitive' first step
  // which is likely to assist in more performant 'stepping'.

  G4double stepTrial;

  //#if 1

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

  //#else

  //  if ( dChordStep > 1000. * fDeltaChord )
  //  {
  //        stepTrial= stepTrialOld * 0.03;   
  //  }
  //  else
  //  {
  //     if ( dChordStep > 100. * fDeltaChord )
  //     {
  //        stepTrial= stepTrialOld * 0.1;   
  //     }
  //     else  // Keep halving the length until dChordStep OK
  //     {
  //        stepTrial= stepTrialOld * 0.5;   
  //     }
  //  }
  //
  //#endif 

  // A more sophisticated chord-finder could figure out a better
  // stepTrial, from dChordStep and the required d_geometry
  //   e.g.
  //      Calculate R, r_helix (eg at orig point)
  //      if( stepTrial < 2 pi  R )
  //          stepTrial = R arc_cos( 1 - fDeltaChord / r_helix )
  //      else    
  //          ??

  return stepTrial;
}


// ...........................................................................

FQUALIFIER
GPFieldTrack
GPChordFinder_ApproxCurvePointS( GPChordFinder *This, 
				 const GPFieldTrack&  CurveA_PointVelocity, 
				 const GPFieldTrack&  CurveB_PointVelocity, 
				 const GPFieldTrack&  ApproxCurveV,
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

  //  GPFieldTrack EndPoint(CurveA_PointVelocity);
  GPFieldTrack EndPoint = CurveA_PointVelocity;
  if(!first){EndPoint= ApproxCurveV;}

  GPThreeVector Point_A,Point_B;
  Point_A=GPFieldTrack_GetPosition(&CurveA_PointVelocity);
  Point_B=GPFieldTrack_GetPosition(&CurveB_PointVelocity);

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
	GPChordFinder_ApproxCurvePointV(This, 
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
    G4double test_step = GPChordFinder_InvParabolic(xa,ya,xb,yb,xc,yc);
    G4double curve;
    if(first)
    {
      curve=abs(GPFieldTrack_GetCurveLength(&EndPoint)
                    -GPFieldTrack_GetCurveLength(&ApproxCurveV));
    }
    else
    {
      test_step=(test_step-xb);
      curve=abs(GPFieldTrack_GetCurveLength(&EndPoint)
                    -GPFieldTrack_GetCurveLength(&CurveB_PointVelocity));
      xb=GPThreeVector_mag(GPThreeVector_sub(CurrentF_Point,Point_B));
    }
      
    if(test_step<=0)    { test_step=0.1*xb; }
    if(test_step>=xb)   { test_step=0.5*xb; }
    if(test_step>=curve){ test_step=0.5*curve; } 

    if(curve*(1.+eps_step)<xb) // Similar to ReEstimate Step from
    {                          // G4VIntersectionLocator
      test_step=0.5*curve;
    }

    GPMagInt_Driver_AccurateAdvance(This->fIntgrDriver,
				    EndPoint,test_step, eps_step, 0);
      
  }
  return EndPoint;
}


// ...........................................................................

FQUALIFIER
GPFieldTrack 
GPChordFinder_ApproxCurvePointV( GPChordFinder *This, 
				 const GPFieldTrack& CurveA_PointVelocity, 
				 const GPFieldTrack& CurveB_PointVelocity, 
				 const GPThreeVector& CurrentE_Point,
				 G4double eps_step)
{
  // If r=|AE|/|AB|, and s=true path lenght (AB)
  // return the point that is r*s along the curve!
 
  GPFieldTrack   Current_PointVelocity = CurveA_PointVelocity; 

  GPThreeVector  CurveA_Point= GPFieldTrack_GetPosition(&CurveA_PointVelocity);
  GPThreeVector  CurveB_Point= GPFieldTrack_GetPosition(&CurveB_PointVelocity);

  GPThreeVector  ChordAB_Vector= GPThreeVector_sub(CurveB_Point,CurveA_Point);
  GPThreeVector  ChordAE_Vector= GPThreeVector_sub(CurrentE_Point,CurveA_Point);

  G4double       ABdist= GPThreeVector_mag(ChordAB_Vector);
  G4double  curve_length;  //  A curve length  of AB
  G4double  AE_fraction; 
  
  curve_length= GPFieldTrack_GetCurveLength(&CurveB_PointVelocity)
              - GPFieldTrack_GetCurveLength(&CurveA_PointVelocity);  
 
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
    GPMagInt_Driver_AccurateAdvance(This->fIntgrDriver,
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


// ......................................................................

FQUALIFIER
void GPChordFinder_PrintStatistics( GPChordFinder *This )
{
  // Print Statistics

  //  G4cout << "GPChordFinder statistics report: " << G4endl;
  //  G4cout 
  //    << "  No trials: " << fTotalNoTrials_FNC
  //    << "  No Calls: "  << fNoCalls_FNC
  //    << "  Max-trial: " <<  fmaxTrials_FNC
  //    << G4endl; 
  //  G4cout 
  //    << "  Parameters: " 
  //    << "  fFirstFraction "  << fFirstFraction
  //    << "  fFractionLast "   << fFractionLast
  //    << "  fFractionNextEstimate " << fFractionNextEstimate
  //    << G4endl; 
}


// ...........................................................................

FQUALIFIER
void GPChordFinder_TestChordPrint( //GPChordFinder *This,
				   G4int    noTrials, 
				   G4int    lastStepTrial, 
				   G4double dChordStep, 
				   G4double nextStepTrial )
{
  ;
  //     G4int oldprec= G4cout.precision(5);
  //     G4cout << " ChF/fnc: notrial " << std::setw( 3) << noTrials 
  //            << " this_step= "       << std::setw(10) << lastStepTrial;
  //  if( std::fabs( (dChordStep / This->fDeltaChord) - 1.0 ) < 0.001 )
  //    {
  //      G4cout.precision(8);
  //    }
  //  else
  //    {
  //      G4cout.precision(6);
  //   }
  //  G4cout << " dChordStep=  " << std::setw(12) << dChordStep;
  //  if( dChordStep > fDeltaChord ) { G4cout << " d+"; }
  //  else                           { G4cout << " d-"; }
  //  G4cout.precision(5);
  //  G4cout <<  " new_step= "       << std::setw(10)
  //	 << fLastStepEstimate_Unconstrained
  //	 << " new_step_constr= " << std::setw(10)
  //	 << lastStepTrial << G4endl;
  //  G4cout << " nextStepTrial = " << std::setw(10) << nextStepTrial << G4endl;
  //  G4cout.precision(oldprec);
}

FQUALIFIER 
void GPChordFinder_SetIntegrationDriver(GPChordFinder *This,
					GPMagInt_Driver* IntegrationDriver)
{
  This->fIntgrDriver=IntegrationDriver;
}

FQUALIFIER
GPMagInt_Driver* GPChordFinder_GetIntegrationDriver(GPChordFinder *This)
{
  return This->fIntgrDriver;
}

FQUALIFIER
G4bool GPChordFinder_AcceptableMissDist(GPChordFinder *This,
					G4double dChordStep)
{ 
  return (dChordStep <= This->fDeltaChord) ;
}

FQUALIFIER
void GPChordFinder_SetChargeMomentumMass(GPChordFinder *This,
					 G4double pCharge, // in e+ units
					 G4double pMomentum,
					 G4double pMass)
{
  GPMagInt_Driver_SetChargeMomentumMass(This->fIntgrDriver,
					pCharge, pMomentum, pMass);
}

FQUALIFIER
G4double GPChordFinder_GetDeltaChord(GPChordFinder *This)
{
  return This->fDeltaChord;
}

FQUALIFIER
void GPChordFinder_SetDeltaChord(GPChordFinder *This,G4double newval)
{
  This->fDeltaChord=newval;
}
// ......................................................................

FQUALIFIER 
G4double GPChordFinder_GetLastStepEstimateUnc(GPChordFinder *This)
{
  return This->fLastStepEstimate_Unconstrained;   
} 

FQUALIFIER 
void GPChordFinder_SetLastStepEstimateUnc(GPChordFinder *This,
					  G4double stepEst )
{
  This->fLastStepEstimate_Unconstrained = stepEst;    
}

FQUALIFIER
void GPChordFinder_ResetStepEstimate(GPChordFinder *This)
{
  This->fLastStepEstimate_Unconstrained = DBL_MAX;    
}

// ......................................................................
FQUALIFIER 
G4int GPChordFinder_GetNoCalls(GPChordFinder *This)     
{ 
  return This->fNoCalls_FNC; 
}

FQUALIFIER 
G4int GPChordFinder_GetNoTrials(GPChordFinder *This)    
{ return This->fTotalNoTrials_FNC; }

FQUALIFIER 
G4int GPChordFinder_GetNoMaxTrials(GPChordFinder *This) 
{ return This->fmaxTrials_FNC; } 

FQUALIFIER 
G4double GPChordFinder_GetFirstFraction(GPChordFinder *This) 
{ return This->fFirstFraction; } 

FQUALIFIER 
G4double GPChordFinder_GetFractionLast(GPChordFinder *This)  
{ return This->fFractionLast; } 

FQUALIFIER 
G4double GPChordFinder_GetFractionNextEstimate(GPChordFinder *This) 
{ return This->fFractionNextEstimate; } 

FQUALIFIER 
G4double GPChordFinder_GetMultipleRadius(GPChordFinder *This)
{ return This->fMultipleRadius; } 

FQUALIFIER 
void GPChordFinder_SetFirstFraction(GPChordFinder *This,
				    G4double val)
{ This->fFirstFraction=val; }

FQUALIFIER 
G4int GPChordFinder_SetVerbose(GPChordFinder *This, G4int newvalue )
{ 
  G4int oldval= This->fStatsVerbose; 
  This->fStatsVerbose = newvalue;
  return oldval; 
}

FQUALIFIER  
void GPChordFinder_AccumulateStatistics(GPChordFinder *This, 
					G4int noTrials ) 
{
  // Statistics 
  This->fTotalNoTrials_FNC += noTrials; 
  This->fNoCalls_FNC++; 
  // if( noTrials >= fmaxTrials_FNC ){
  if (noTrials > This->fmaxTrials_FNC ) { 
    This->fmaxTrials_FNC=noTrials; 
    // fnoTimesMaxTrFNC=0; 
  } else { 
    // fnoTimesMaxTrFNC++; 
  } 
  // } 
}

// A  member that calculates the inverse parabolic through
// the three points (x,y) and returns the value x that, for the
// inverse parabolic, corresponds to y=0.
FQUALIFIER 
G4double GPChordFinder_InvParabolic ( const G4double xa, 
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
