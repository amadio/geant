
//
// $Id: GVIntegratorDriver.cxx
//
// Implementation for class GUIntegrationDriver
//   Tracking in space dependent magnetic field
//
// History of major changes:
//  8 Dec 14  J. Apostolakis:   First version
// --------------------------------------------------------------------

#include <iostream>
#include <iomanip>

// #include "globals.hh"
// #include "G4SystemOfUnits.hh"
// #include "G4GeometryTolerance.hh"
#include "GUFieldTrack.h"
#include "GUIntegrationDriver.h"

//  Stepsize can increase by no more than 5.0
//           and decrease by no more than 1/10. = 0.1
//
const double GUIntegrationDriver::max_stepping_increase = 5.0;
const double GUIntegrationDriver::max_stepping_decrease = 0.1;

const double perMillion  = 1.0e-6; 
const double perThousand = 1.0e-3;

//  The (default) maximum number of steps is Base
//  divided by the order of Stepper
//
const int  GUIntegrationDriver::fMaxStepBase = 250;  // Was 5000

#ifndef G4NO_FIELD_STATISTICS
#define GVFLD_STATS  1
#endif

// To add much printing for debugging purposes, uncomment the following
// and set verbose level to 1 or higher value !
// #define  GUDEBUG_FIELD 1    

// ---------------------------------------------------------

//  Constructor
//
GUIntegrationDriver::GUIntegrationDriver( double                hminimum, 
                                  GUVIntegrationStepper *pStepper,
                                  int                   numComponents,
                                  int                   statisticsVerbose)
   : fMinimumStep( hminimum ),
     fSmallestFraction( 1.0e-12 ), 
     fNoIntegrationVariables(numComponents),
     fMinNoVars(12),
     fNoVars( std::max( fNoIntegrationVariables, fMinNoVars )),
     fSafetyFactor(0.9),
     fPowerShrink(0.0),   //  exponent for shrinking
     fPowerGrow(0.0),    //  exponent for growth
     fErrcon(0.0),
     fSurfaceTolerance( 1.0e-6),
     fStatisticsVerboseLevel(statisticsVerbose),
     fNoTotalSteps(0),  fNoBadSteps(0), fNoSmallSteps(0),
     fNoInitialSmallSteps(0), 
     fDyerr_max(0.0), fDyerr_mx2(0.0), 
     fDyerrPos_smTot(0.0), fDyerrPos_lgTot(0.0), fDyerrVel_lgTot(0.0), 
     fSumH_sm(0.0), fSumH_lg(0.0),
     fVerboseLevel(0)
{  
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12

  RenewStepperAndAdjust( pStepper );
  fMaxNoSteps = fMaxStepBase / fpStepper->IntegratorOrder();
#ifdef GUDEBUG_FIELD
  fVerboseLevel=2;
#endif

  if( (fVerboseLevel > 0) || (fStatisticsVerboseLevel > 1) )
  {
     std::cout << "MagIntDriver version: Accur-Adv: "
           << "invE_nS, QuickAdv-2sqrt with Statistics "
#ifdef GVFLD_STATS
           << " enabled "
#else
           << " disabled "
#endif
           << std::endl;
  }
}

//  Copy Constructor - used by Clone
//
GUIntegrationDriver::GUIntegrationDriver( const GUIntegrationDriver& right ) 
   : fMinimumStep( right.fMinimumStep ),
     fSmallestFraction( right.fSmallestFraction ),
     fNoIntegrationVariables( right.fNoIntegrationVariables ),
     fMinNoVars( right.fMinNoVars ),
     fNoVars( std::max( fNoIntegrationVariables, fMinNoVars )),
     fSafetyFactor( right.fSafetyFactor ),
     fPowerShrink( right.fPowerShrink ),
     fPowerGrow(   right.fPowerGrow),
     fErrcon( right.fErrcon ),
     fSurfaceTolerance( right.fSurfaceTolerance ),
     fStatisticsVerboseLevel( right.fStatisticsVerboseLevel ),
     fNoTotalSteps(0),  fNoBadSteps(0), fNoSmallSteps(0),
     fNoInitialSmallSteps(0),
     fDyerr_max(0.0), fDyerr_mx2(0.0), 
     fDyerrPos_smTot(0.0), fDyerrPos_lgTot(0.0), fDyerrVel_lgTot(0.0), 
     fSumH_sm(0.0), fSumH_lg(0.0),
     fVerboseLevel( right.fVerboseLevel )
{  
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
   const GUVIntegrationStepper *protStepper = right.GetStepper();
   fpStepper= protStepper->Clone();
   
   RenewStepperAndAdjust( fpStepper );
   fMaxNoSteps = fMaxStepBase / fpStepper->IntegratorOrder();

  if( (fVerboseLevel > 0) || (fStatisticsVerboseLevel > 1) )
  {
     std::cout << "MagIntDriver version: Accur-Adv: "
           << "invE_nS, QuickAdv-2sqrt with Statistics "
#ifdef GVFLD_STATS
           << " enabled "
#else
           << " disabled "
#endif
           << std::endl;
  }
}


// ---------------------------------------------------------

//  Destructor
//
GUIntegrationDriver::~GUIntegrationDriver()
{ 
  if( fStatisticsVerboseLevel > 1 )
  {
    PrintStatisticsReport();
  }
}


GUIntegrationDriver* GUIntegrationDriver::Clone() const
{
   return new GUIntegrationDriver(*this);
}

// ---------------------------------------------------------

bool
GUIntegrationDriver::AccurateAdvance(const GUFieldTrack& yInput,
                                              double     hstep,
                                              double     eps,
                                           GUFieldTrack& yOutput,
                                              double     hinitial
   )
{
  // Runge-Kutta driver with adaptive stepsize control. Integrate starting
  // values at y_current over hstep x2 with accuracy eps. 
  // On output ystart is replaced by values at the end of the integration 
  // interval. RightHandSide is the right-hand side of ODE system. 
  // The source is similar to odeint routine from NRC p.721-722 .

  int nstp, i, no_warnings=0;
  double x, hnext, hdid, h;
  // double charge= y_current.GetCharge(); 

#ifdef GUDEBUG_FIELD
  static int dbg=1;
  static int nStpPr=50;   // For debug printing of long integrations
  double ySubStepStart[GUFieldTrack::ncompSVEC];
  // GUFieldTrack  yFldTrkStart(y_current);
#endif

  double y[GUFieldTrack::ncompSVEC], dydx[GUFieldTrack::ncompSVEC];
  double ystart[GUFieldTrack::ncompSVEC], yEnd[GUFieldTrack::ncompSVEC]; 
  double  x1, x2;
  bool succeeded = true, lastStepSucceeded;

  double startCurveLength;

  int  noFullIntegr=0, noSmallIntegr = 0 ;
  // G4ThreadLocal
  static int  noGoodSteps =0 ;  // Bad = chord > curve-len 
  const  int  nvar= fNoVars;

  // GUFieldTrack yStartFT(yInput);

  //  Ensure that hstep > 0
  //
  if( hstep <= 0.0 )
  { 
    if(hstep==0.0)
    {
      std::cerr << "Proposed step is zero; hstep = " << hstep << " !";
      return succeeded; 
    }
    else
    { 
      std::cerr << "Invalid run condition." << std::endl
                << "Proposed step is negative; hstep = " << hstep << "." << std::endl;
      return false;
    }
  }

  yInput.DumpToArray( ystart );

  startCurveLength= yInput.GetCurveLength();
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

  bool lastStep= false;
  nstp=1;

  // fpStepper->InitializeCharge( charge );   // May be added
  // OLD: fpStepper->GetEquationOfMotion()->InitializeCharge( charge );
  
  do
  {
    ThreeVector StartPos( y[0], y[1], y[2] );

#ifdef GUDEBUG_FIELD
    double xSubStepStart= x; 
    for (i=0;i<nvar;i++)  { ySubStepStart[i] = y[i]; }
    // yFldTrkStart.LoadFromArray(y, fNoIntegrationVariables);
    // yFldTrkStart.SetCurveLength(x);
#endif

    // Old method - inline call to Equation of Motion
    //   fpStepper->RightHandSide( y, dydx );
    // New method allows to cache field, or state (eg momentum magnitude)
    // fpStepper->ComputeRightHandSide( y, charge, dydx );

    // Back to simple, old method   - JA. 16 Oct 2015
    fpStepper->RightHandSideVIS( y, /*charge,*/ dydx );   // TODO: change to inline
    fNoTotalSteps++;

    // Perform the Integration
    //      
    if( h > fMinimumStep )
    { 
      OneGoodStep(y,dydx,x,h,eps,hdid,hnext) ;
      //--------------------------------------
      lastStepSucceeded= (hdid == h);   
#ifdef GUDEBUG_FIELD
      if (dbg>2) {
        PrintStatus( ySubStepStart, xSubStepStart, y, x, h,  nstp); // Only
      }
#endif
    }
    else
    {
      GUFieldTrack yFldTrk( ThreeVector(0,0,0), 
                            ThreeVector(0,0,0), 0., 0. );
      double dchord_step, dyerr_len_sq, dyerr_mom_rel_sq;   // What to do with these ?
      yFldTrk.LoadFromArray(y, fNoIntegrationVariables); 
      yFldTrk.SetCurveLength( x );

      QuickAdvance( yFldTrk, dydx, h, dchord_step, dyerr_len_sq, dyerr_mom_rel_sq ); 
      //----------------------------------------------------------------

      yFldTrk.DumpToArray(y);    


#ifdef GVFLD_STATS
      fNoSmallSteps++;
      double dyerr_len = std::sqrt(dyerr_len_sq); 
      if ( dyerr_len > fDyerr_max)  { fDyerr_max= dyerr_len; }
      fDyerrPos_smTot += dyerr_len;
      fSumH_sm += h;  // Length total for 'small' steps
      if (nstp<=1)  { fNoInitialSmallSteps++; }
#endif
#ifdef GUDEBUG_FIELD
      if (dbg>1)
      {
        if(fNoSmallSteps<2) { PrintStatus(ySubStepStart, x1, y, x, h, -nstp); }
        std::cout << "Another sub-min step, no " << fNoSmallSteps 
               << " of " << fNoTotalSteps << " this time " << nstp << std::endl; 
        PrintStatus( ySubStepStart, x1, y, x, h,  nstp);   // Only this
        std::cout << " dyerr= " << dyerr_len << " relative = " << dyerr_len / h 
               << " epsilon= " << eps << " hstep= " << hstep 
               << " h= " << h << " hmin= " << fMinimumStep << std::endl;
      }
#endif        
      if( h == 0.0 ){
         std::cerr << "ERROR in G4UIntegationDriver::AccurateAdvance: integration failure. " << std::endl;
         std::cerr << "Integration Step became Zero!" << std::endl;
         exit(1); 
      }
      double dyerr_sq = std::min( dyerr_len_sq / (h*h) , dyerr_mom_rel_sq ); 
      double dyerr = std::sqrt(dyerr_sq); 
      // dyerr = dyerr_len / h;
      hdid= h;
      x += hdid;

      // Compute suggested new step
      hnext= ComputeNewStepSize( dyerr/eps, h);

      // .. hnext= ComputeNewStepSize_WithinLimits( dyerr/eps, h);
      lastStepSucceeded= (dyerr<= eps);
    }

    if (lastStepSucceeded)  { noFullIntegr++; }
    else                    { noSmallIntegr++; }

    ThreeVector EndPos( y[0], y[1], y[2] );

#ifdef  GUDEBUG_FIELD
    if( (dbg>0) && (dbg<=2) && (nstp>nStpPr))
    {
      if( nstp==nStpPr )  { std::cout << "***** Many steps ****" << std::endl; }
      std::cout << "MagIntDrv: " ; 
      std::cout << "hdid="  << std::setw(12) << hdid  << " "
             << "hnext=" << std::setw(12) << hnext << " " 
	     << "hstep=" << std::setw(12) << hstep << " (requested) " 
	     << std::endl;
      PrintStatus( ystart, x1, y, x, h, (nstp==nStpPr) ? -nstp: nstp); 
    }
#endif

    // Check the endpoint
    double endPointDist= (EndPos-StartPos).Mag(); 
    if ( endPointDist >= hdid*(1.+perMillion) )
    {
      fNoBadSteps++;

      // Issue a warning only for gross differences -
      // we understand how small difference occur.
      if ( endPointDist >= hdid*(1.+perThousand) )
      { 
#ifdef GUDEBUG_FIELD
        if (dbg)
        {
          WarnEndPointTooFar ( endPointDist, hdid, eps, dbg ); 
          std::cerr << "  Total steps:  bad " << fNoBadSteps
                 << " good " << noGoodSteps << " current h= " << hdid
                 << std::endl;
          PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);  
        }
#endif
        no_warnings++;
      }
    }
    else
    {
      noGoodSteps ++;
    } 
// #endif

    //  Avoid numerous small last steps
    if( (h < eps * hstep) || (h < fSmallestFraction * startCurveLength) )
    {
      // No more integration -- the next step will not happen
      lastStep = true;  
    }
    else
    {
      // Check the proposed next stepsize
      if(std::fabs(hnext) <= fMinimumStep)
      {
#ifdef  GUDEBUG_FIELD
        // If simply a very small interval is being integrated, do not warn
        if( (x < x2 * (1-eps) ) &&        // The last step can be small: OK
            (std::fabs(hstep) > Hmin()) ) // and if we are asked, it's OK
        {
          if(dbg>0)
          {
            WarnSmallStepSize( hnext, hstep, h, x-x1, nstp );  
            PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);
          }
          no_warnings++;
        }
#endif
        // Make sure that the next step is at least Hmin.
        h = fMinimumStep;
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
#ifdef GUDEBUG_FIELD
        if (dbg>2)
        {
          int prec= std::cout.precision(12); 
          std::cout << "Warning: GVIntegratorDriver::AccurateAdvance"
                 << std::endl
                 << "  Integration step 'h' became "
                 << h << " due to roundoff. " << std::endl
		 << " Calculated as difference of x2= "<< x2 << " and x=" << x
                 << "  Forcing termination of advance." << std::endl;
          std::cout.precision(prec);
        }          
#endif
      }
    }
  } while ( ((nstp++)<=fMaxNoSteps) && (x < x2) && (!lastStep) );
     // Have we reached the end ?
     // --> a better test might be x-x2 > an_epsilon

  succeeded=  (x>=x2);  // If it was a "forced" last step

  for (i=0;i<nvar;i++)  { yEnd[i] = y[i]; }

  // Put back the values.
  yOutput.LoadFromArray( yEnd, fNoIntegrationVariables );
  yOutput.SetCurveLength( x );

  if(nstp > fMaxNoSteps)
  {
    no_warnings++;
    succeeded = false;
#ifdef GUDEBUG_FIELD
    if (dbg)
    {
      WarnTooManyStep( x1, x2, x );  //  Issue WARNING
      PrintStatus( yEnd, x1, y, x, hstep, -nstp);
    }
#endif
  }

#ifdef GUDEBUG_FIELD
  if( dbg && no_warnings )
  {
    std::cerr << "GVIntegratorDriver exit status: no-steps " << nstp <<std::endl;
    PrintStatus( yEnd, x1, y, x, hstep, nstp);
  }
#endif

  // No longer relevant -- IntDriver does not control charge anymore
  // OLD: Finished for now - value of charge is 'revoked'
  //  fpStepper->GetEquationOfMotion()->InformDone();  
  //  fpStepper->InformDone();  // May be added
  
  return succeeded;
}  // end of AccurateAdvance ...........................

// ---------------------------------------------------------

void
GUIntegrationDriver::WarnSmallStepSize( double hnext, double hstep, 
                                    double h, double xDone,
                                    int nstp)
{
  static int noWarningsIssued =0;   // thread_local
  const  int maxNoWarnings =  10;   // Number of verbose warnings
  // std::ostringstream message;
  // typedef std::cerr message;
  std::cerr << " WARNING from GUIntegrationDriver::WarnSmallStepSize() " << std::endl;
  if( (noWarningsIssued < maxNoWarnings) || fVerboseLevel > 10 )
  {
    std::cerr << "The stepsize for the next iteration, " << hnext
            << ", is too small - in Step number " << nstp << "." << std::endl
            << "The minimum for the driver is " << GetHmin()  << std::endl
            << "Requested integr. length was " << hstep << " ." << std::endl
            << "The size of this sub-step was " << h     << " ." << std::endl
            << "The integrations has already gone " << xDone;
  }
  else
  {
    std::cerr << "Too small 'next' step " << hnext
            << ", step-no: " << nstp << std::endl
            << ", this sub-step: " << h     
            << ",  req_tot_len: " << hstep 
            << ", done: " << xDone << ", min: " << GetHmin();
  }
  // G4Exception("GUIntegrationDriver::WarnSmallStepSize()", "GeomField1001",
  //             JustWarning, message);
  noWarningsIssued++;
}

// ---------------------------------------------------------

void
GUIntegrationDriver::WarnTooManyStep( double x1start, 
                                  double x2end, 
                                  double xCurrent)
{
   // std::ostringstream message;
   std::cerr << "WARNING from GUIntegrationDriver::WarnTooManyStep()" << std::endl;
   std::cerr << "The number of steps used in the Integration driver"
             << " (Runge-Kutta) is too many." << std::endl
             << "Integration of the interval was not completed !" << std::endl
             << "Only a " << (xCurrent-x1start)*100/(x2end-x1start)
             << " % fraction of it was done." << std::endl;
   // G4Exception("GUIntegrationDriver::WarnTooManyStep()", "GeomField1001",
   //             JustWarning, message);
}

// ---------------------------------------------------------

void
GUIntegrationDriver::WarnEndPointTooFar (double endPointDist, 
                                     double   h , 
                                     double  eps,
                                     int     dbg)
{
  static  double maxRelError=0.0; // thread_local
  bool isNewMax, prNewMax;

  isNewMax = endPointDist > (1.0 + maxRelError) * h;
  prNewMax = endPointDist > (1.0 + 1.05 * maxRelError) * h;
  if( isNewMax ) { maxRelError= endPointDist / h - 1.0; }

  if( dbg && (h > fSurfaceTolerance) 
          && ( (dbg>1) || prNewMax || (endPointDist >= h*(1.+eps) ) ) )
  { 
    static int noWarnings = 0;  // thread_local
    // std::ostringstream message;
    std::cerr << "WARNING in GUIntegrationDriver::WarnEndPointTooFar()" << std::endl;
    if( (noWarnings ++ < 10) || (dbg>2) )
    {
      std::cerr << "The integration produced an end-point which " << std::endl
              << "is further from the start-point than the curve length."
              << std::endl;
    }
    std::cerr << "  Distance of endpoints = " << endPointDist
              << ", curve length = " << h << std::endl
              << "  Difference (curveLen-endpDist)= " << (h - endPointDist)
              << ", relative = " << (h-endPointDist) / h 
              << ", epsilon =  " << eps << std::endl;
    // G4Exception("GUIntegrationDriver::WarnEndPointTooFar()", "GeomField1001",
    //             JustWarning, message);
  }
}

// ---------------------------------------------------------

void
GUIntegrationDriver::OneGoodStep(      double y[],        // InOut
                             const double dydx[],
                                   double& x,         // InOut
                                   double htry,
                                   double eps_rel_max,
                                   double& hdid,      // Out
                                   double& hnext )    // Out

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
  double errmax_sq;
  double h, htemp, xnew ;

  double yerr[GUFieldTrack::ncompSVEC], ytemp[GUFieldTrack::ncompSVEC];

  h = htry ; // Set stepsize to the initial trial value

  double inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  double errpos_sq=0.0;    // square of displacement error
  double errmom_sq=0.0;    // square of momentum vector difference
  //double errspin_sq=0.0;   // square of spin vector difference

  int iter;

  static int tot_no_trials=0;  // thread_local
  const int max_trials=100; 

  // ThreeVector Spin(y[9],y[10],y[11]);
  // double   spin_mag2 =Spin.Mag2() ;
  // bool     hasSpin= (spin_mag2 > 0.0); 

  for (iter=0; iter<max_trials ;iter++)
  {
    tot_no_trials++;
    fpStepper-> StepWithErrorEstimate(y,dydx,h,ytemp,yerr);
    //          *********************
    double eps_pos = eps_rel_max * std::max(h, fMinimumStep);  // Uses remaining step 'h'
    double inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos);

    // Evaluate accuracy
    //
    errpos_sq =  yerr[0]*yerr[0] + yerr[1]*yerr[1] + yerr[2]*yerr[2] ; 
    errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

    // Accuracy for momentum
    double magmom_sq=  y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    double sumerr_sq=  yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5]; 
    if( magmom_sq > 0.0 ) { 
       // double inv_magmom_sq = 1.0 / magmom_sq; 
       // errmom_sq = sumerr_sq * inv_magmom_sq;
       errmom_sq = sumerr_sq / magmom_sq; 
    } else {
       std::cerr << "** GUIntegrationDriver: found case of zero momentum." 
                 << " iteration=  " << iter << " h= " << h << std::endl; 
       errmom_sq = sumerr_sq;
    }
    errmom_sq *= inv_eps_vel_sq;
    errmax_sq = std::max( errpos_sq, errmom_sq ); // Square of maximum error

    // if( hasSpin )
    // { 
    //    // Accuracy for spin
    //     errspin_sq =  ( sqr(yerr[9]) + sqr(yerr[10]) + sqr(yerr[11]) )
    //                     /  spin_mag2; // ( sqr(y[9]) + sqr(y[10]) + sqr(y[11]) );
    //   errspin_sq *= inv_eps_vel_sq;
    //   errmax_sq = std::max( errmax_sq, errspin_sq ); 
    // }

    if ( errmax_sq <= 1.0 )  { break; } // Step succeeded. 

    // Step failed; compute the size of retrial Step.
    htemp = fSafetyFactor *h* std::pow( errmax_sq, 0.5*fPowerShrink );

    if (htemp >= 0.1*h)  { h = htemp; }  // Truncation error too large,
    else  { h = 0.1*h; }                 // reduce stepsize, but no more
                                         // than a factor of 10
    xnew = x + h;
    if(xnew == x)
    {
      std::cerr << "GVIntegratorDriver::OneGoodStep:" << std::endl
             << "  Stepsize underflow in Stepper " << std::endl ;
      std::cerr << "  Step's start x=" << x << " and end x= " << xnew 
             << " are equal !! " << std::endl
             <<"  Due to step-size= " << h 
             << " . Note that input step was " << htry << std::endl;
      break;
    }
  }

#ifdef GVFLD_STATS
  // Sum of squares of position error // and momentum dir (underestimated)
  fSumH_lg += h; 
  fDyerrPos_lgTot += errpos_sq;
  fDyerrVel_lgTot += errmom_sq * h * h; 
#endif

  // Compute size of next Step
  if (errmax_sq > fErrcon*fErrcon)
  { 
    hnext = GetSafety()*h*std::pow(errmax_sq, 0.5*GetPowerGrow());
  }
  else
  {
    hnext = max_stepping_increase*h ; // No more than a factor of 5 increase
  }
  x += (hdid = h);

  for(int k=0;k<fNoIntegrationVariables;k++) { y[k] = ytemp[k]; }

  return;
}   // end of  OneGoodStep .............................

//----------------------------------------------------------------------

#define SQR(a)   ((a)*(a))

// QuickAdvance just tries one Step - it does not ensure accuracy
//
bool  GUIntegrationDriver::QuickAdvance(       
                            GUFieldTrack& y_posvel,         // INOUT
                            const double     dydx[],  
                                  double     hstep,       // In
                                  double&    dchord_step,
                                  double&    dyerr_pos_sq,
                                  double&    dyerr_mom_rel_sq )  
{
  // double dyerr_pos_sq, dyerr_mom_rel_sq;  
  double yerr_vec[GUFieldTrack::ncompSVEC],
           yarrin[GUFieldTrack::ncompSVEC], yarrout[GUFieldTrack::ncompSVEC]; 
  double s_start;
  double dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

  static int no_call=0;  // thread_local 
  no_call ++; 

  // Move data into array
  y_posvel.DumpToArray( yarrin );      //  yarrin  <== y_posvel 
  s_start = y_posvel.GetCurveLength();

  // Do an Integration Step
  fpStepper-> StepWithErrorEstimate(yarrin, dydx, hstep, yarrout, yerr_vec) ; 
  //          *********************

  // Estimate curve-chord distance
  dchord_step= fpStepper-> DistChord();
  //                       *********

  // Put back the values.  yarrout ==> y_posvel
  y_posvel.LoadFromArray( yarrout, fNoIntegrationVariables );
  y_posvel.SetCurveLength( s_start + hstep );

#ifdef  GUDEBUG_FIELD
  if(fVerboseLevel>2)
  {
    std::cout << "G4MagIntDrv: Quick Advance" << std::endl;
    PrintStatus( yarrin, s_start, yarrout, s_start+hstep, hstep,  1); 
  }
#endif

  // A single measure of the error   
  //      TO-DO :  account for  energy,  spin, ... ? 
  vel_mag_sq   = ( SQR(yarrout[3])+SQR(yarrout[4])+SQR(yarrout[5]) );
  inv_vel_mag_sq = 1.0 / vel_mag_sq; 
  dyerr_pos_sq = ( SQR(yerr_vec[0])+SQR(yerr_vec[1])+SQR(yerr_vec[2]));
  dyerr_mom_sq = ( SQR(yerr_vec[3])+SQR(yerr_vec[4])+SQR(yerr_vec[5]));
  dyerr_mom_rel_sq = dyerr_mom_sq * inv_vel_mag_sq;

  // Calculate also the change in the momentum squared also ???
  // double veloc_square = y_posvel.GetVelocity().mag2();
  // ...

#ifdef RETURN_A_NEW_STEP_LENGTH
  // The following step cannot be done here because "eps" is not known.
  dyerr_len = std::sqrt( dyerr_len_sq ); 
  dyerr_len_sq /= eps ;

  // Look at the velocity deviation ?
  //  SQR(yerr_vec[3])+SQR(yerr_vec[4])+SQR(yerr_vec[5]));

  // Set suggested new step
  hstep= ComputeNewStepSize( dyerr_len, hstep);
#endif

#if 0
  if( dyerr_pos_sq > ( dyerr_mom_rel_sq * SQR(hstep) ) )
  {
    dyerr = std::sqrt(dyerr_pos_sq);
  }
  else
  {
    // Scale it to the current step size - for now
    dyerr = std::sqrt(dyerr_mom_rel_sq) * hstep;
  }
#endif

  return true;
}

// --------------------------------------------------------------------------

#ifdef QUICK_ADV_ARRAY_IN_AND_OUT
bool  GUIntegrationDriver::QuickAdvance(       
                                  double     yarrin[],    // In
                            const double     dydx[],  
                                  double     hstep,       // In
                                  double     yarrout[],
                                  double&    dchord_step,
                                  double&    dyerr )      // In length
{
   std::cerr << "ERROR in GUIntegrationDriver::QuickAdvance()" << std::endl;
   std::cerr << "      Method is not yet implemented." << std::endl;

   //            FatalException, "Not yet implemented.");
   dyerr = dchord_step = hstep * yarrin[0] * dydx[0];
   yarrout[0]= yarrin[0];
   exit(1); 
}
#endif 

// --------------------------------------------------------------------------

//  This method computes new step sizes - but does not limit changes to
//   within  certain factors
// 
double 
GUIntegrationDriver::ComputeNewStepSize( 
                          double  errMaxNorm,    // max error  (normalised)
                          double  hstepCurrent)  // current step size
{
  double hnew;

  // Compute size of next Step for a failed step
  if(errMaxNorm > 1.0 )
  {
    // Step failed; compute the size of retrial Step.
    hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPowerShrink()) ;
  } else if(errMaxNorm > 0.0 ) {
    // Compute size of next Step for a successful step
    hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPowerGrow()) ;
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
double 
GUIntegrationDriver::ComputeNewStepSize_WithinLimits( 
                          double  errMaxNorm,    // max error  (normalised)
                          double  hstepCurrent)  // current step size
{
  double hnew;

  // Compute size of next Step for a failed step
  if (errMaxNorm > 1.0 )
  {
    // Step failed; compute the size of retrial Step.
    hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPowerShrink()) ;
  
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
    if (errMaxNorm > fErrcon)
     { hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPowerGrow()); }
    else  // No more than a factor of 5 increase
     { hnew = max_stepping_increase * hstepCurrent; }
  }
  return hnew;
}

// ---------------------------------------------------------------------------

void GUIntegrationDriver::PrintStatus( const double*   StartArr,  
                                   double          xstart,
                                   const double*   CurrentArr, 
                                   double          xcurrent,
                                   double          requestStep, 
                                   int             subStepNo)
  // Potentially add as arguments:  
  //                                 <dydx>           - as Initial Force
  //                                 stepTaken(hdid)  - last step taken
  //                                 nextStep (hnext) - proposal for size
{
   GUFieldTrack  StartFT(ThreeVector(0,0,0),
                 ThreeVector(0,0,0), 0., 0. );
   GUFieldTrack  CurrentFT (StartFT);

   StartFT.LoadFromArray( StartArr, fNoIntegrationVariables); 
   StartFT.SetCurveLength( xstart);
   CurrentFT.LoadFromArray( CurrentArr, fNoIntegrationVariables); 
   CurrentFT.SetCurveLength( xcurrent );

   PrintStatus(StartFT, CurrentFT, requestStep, subStepNo ); 
}

// ---------------------------------------------------------------------------

void GUIntegrationDriver::PrintStatus(
                  const GUFieldTrack&  StartFT,
                  const GUFieldTrack&  CurrentFT, 
                  double             requestStep, 
                  int                subStepNo)
{
    int verboseLevel= fVerboseLevel;
    static int noPrecision= 5;  // thread_local ? 
    int oldPrec= std::cout.precision(noPrecision);
    // std::cout.setf(ios_base::fixed,ios_base::floatfield);

  // const ThreeVector StartPosition=       StartFT.GetPosition();
    const ThreeVector StartUnitVelocity=   StartFT.GetMomentumDirection();
  // const ThreeVector CurrentPosition=     CurrentFT.GetPosition();
    const ThreeVector CurrentUnitVelocity= CurrentFT.GetMomentumDirection();

    double  DotStartCurrentVeloc= StartUnitVelocity.Dot(CurrentUnitVelocity);

    double step_len= CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
    double subStepSize = step_len;
     
    if( (subStepNo <= 1) || (verboseLevel > 3) )
    {
       subStepNo = - subStepNo;        // To allow printing banner

       std::cout << std::setw( 6)  << " " << std::setw( 25)
              << " GUIntegrationDriver: Current Position  and  Direction" << " "
              << std::endl; 
       std::cout << std::setw( 5) << "Step#" << " "
              << std::setw( 7) << "s-curve" << " "
              << std::setw( 9) << "X(mm)" << " "
              << std::setw( 9) << "Y(mm)" << " "  
              << std::setw( 9) << "Z(mm)" << " "
              << std::setw( 8) << " N_x " << " "
              << std::setw( 8) << " N_y " << " "
              << std::setw( 8) << " N_z " << " "
              << std::setw( 8) << " N^2-1 " << " "
              << std::setw(10) << " N(0).N " << " "
              << std::setw( 7) << "KinEner " << " "
              << std::setw(12) << "Track-l" << " "   // Add the Sub-step ??
              << std::setw(12) << "Step-len" << " " 
              << std::setw(12) << "Step-len" << " " 
              << std::setw( 9) << "ReqStep" << " "  
              << std::endl;
    }

    if( (subStepNo <= 0) )
    {
      PrintStat_Aux( StartFT,  requestStep, 0., 
                       0,        0.0,         1.0);
      //*************
    }

    if( verboseLevel <= 3 )
    {
      std::cout.precision(noPrecision);
      PrintStat_Aux( CurrentFT, requestStep, step_len, 
                     subStepNo, subStepSize, DotStartCurrentVeloc );
      //*************
    }

    else // if( verboseLevel > 3 )
    {
       //  Multi-line output
       
       // std::cout << "Current  Position is " << CurrentPosition << std::endl 
       //    << " and UnitVelocity is " << CurrentUnitVelocity << std::endl;
       // std::cout << "Step taken was " << step_len  
       //    << " out of PhysicalStep= " <<  requestStep << std::endl;
       // std::cout << "Final safety is: " << safety << std::endl;
       // std::cout << "Chord length = " << (CurrentPosition-StartPosition).mag()
       //        << std::endl << std::endl; 
    }
    std::cout.precision(oldPrec);
}

// ---------------------------------------------------------------------------

void GUIntegrationDriver::PrintStat_Aux(
                  const GUFieldTrack&  aGUFieldTrack,
                  double             requestStep, 
                  double             step_len,
                  int                subStepNo,
                  double             subStepSize,
                  double             dotVeloc_StartCurr)
{
    const ThreeVector Position=      aGUFieldTrack.GetPosition();
    const ThreeVector UnitVelocity=  aGUFieldTrack.GetMomentumDirection();
 
    if( subStepNo >= 0)
    {
       std::cout << std::setw( 5) << subStepNo << " ";
    }
    else
    {
       std::cout << std::setw( 5) << "Start" << " ";
    }
    double curveLen= aGUFieldTrack.GetCurveLength();
    std::cout << std::setw( 7) << curveLen;
    std::cout << std::setw( 9) << Position.x() << " "
           << std::setw( 9) << Position.y() << " "
           << std::setw( 9) << Position.z() << " "
           << std::setw( 8) << UnitVelocity.x() << " "
           << std::setw( 8) << UnitVelocity.y() << " "
           << std::setw( 8) << UnitVelocity.z() << " ";
    int oldprec= std::cout.precision(3);
    std::cout << std::setw( 8) << UnitVelocity.Mag2()-1.0 << " ";
    std::cout.precision(6);
    std::cout << std::setw(10) << dotVeloc_StartCurr << " ";
    std::cout.precision(oldprec);
    // std::cout << std::setw( 7) << aGUFieldTrack.GetKineticEnergy();
    std::cout << std::setw(12) << step_len << " ";

    static double oldCurveLength= 0.0;    // thread_local
    static double oldSubStepLength= 0.0; // thread_local
    static int oldSubStepNo= -1;// thread_local

    double subStep_len=0.0;
    if( curveLen > oldCurveLength )
    {
      subStep_len= curveLen - oldCurveLength;
    }
    else if (subStepNo == oldSubStepNo)
    {
      subStep_len= oldSubStepLength;
    }
    oldCurveLength= curveLen;
    oldSubStepLength= subStep_len;

    std::cout << std::setw(12) << subStep_len << " "; 
    std::cout << std::setw(12) << subStepSize << " "; 
    if( requestStep != -1.0 )
    {
      std::cout << std::setw( 9) << requestStep << " ";
    }
    else
    {
       std::cout << std::setw( 9) << " InitialStep " << " ";
    }
    std::cout << std::endl;
}

// ---------------------------------------------------------------------------

void GUIntegrationDriver::PrintStatisticsReport()
{
  int noPrecBig= 6;
  int oldPrec= std::cout.precision(noPrecBig);

  std::cout << "GUIntegrationDriver Statistics of steps undertaken. " << std::endl;
  std::cout << "GUIntegrationDriver: Number of Steps: "
         << " Total= " <<  fNoTotalSteps
         << " Bad= "   <<  fNoBadSteps 
         << " Small= " <<  fNoSmallSteps 
         << " Non-initial small= " << (fNoSmallSteps-fNoInitialSmallSteps)
         << std::endl;

#ifdef GVFLD_STATS
  std::cout << "MID dyerr: " 
         << " maximum= " << fDyerr_max 
         << " Sum small= " << fDyerrPos_smTot 
         << " std::sqrt(Sum large^2): pos= " << std::sqrt(fDyerrPos_lgTot)
         << " vel= " << std::sqrt( fDyerrVel_lgTot )
         << " Total h-distance: small= " << fSumH_sm 
         << " large= " << fSumH_lg
         << std::endl;

#if 0
  int noPrecSmall=4; 
  // Single line precis of statistics ... optional
  std::cout.precision(noPrecSmall);
  std::cout << "MIDnums: " << fMinimumStep
         << "   " << fNoTotalSteps 
         << "  "  <<  fNoSmallSteps
         << "  "  << fNoSmallSteps-fNoInitialSmallSteps
         << "  "  << fNoBadSteps         
         << "   " << fDyerr_max
         << "   " << fDyerr_mx2 
         << "   " << fDyerrPos_smTot 
         << "   " << fSumH_sm
         << "   " << fDyerrPos_lgTot
         << "   " << fDyerrVel_lgTot
         << "   " << fSumH_lg
         << std::endl;
#endif 
#endif 

 std::cout.precision(oldPrec);
}
 
// ---------------------------------------------------------------------------

void GUIntegrationDriver::SetSmallestFraction(double newFraction)
{
  if( (newFraction > 1.e-16) && (newFraction < 1e-8) )
  {
    fSmallestFraction= newFraction;
  }
  else
  { 
    std::cerr << "Warning: SmallestFraction not changed. " << std::endl
           << "  Proposed value was " << newFraction << std::endl
           << "  Value must be between 1.e-8 and 1.e-16" << std::endl;
  }
}
