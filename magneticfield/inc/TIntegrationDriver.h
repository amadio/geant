//
//

#ifndef TIntegrationDriver_Def
#define TIntegrationDriver_Def

#include <iomanip>

#include "GUFieldTrack.hh"
#include "globals.hh"
// #include "SystemOfUnits.hh"
// #include "GeometryTolerance.hh"
#include "TCommonType.hh"

#ifndef NO_FIELD_STATISTICS
#define FIELD_STATS  1
#endif

// To add much printing for debugging purposes, uncomment the following
// and set verbose level to 1 or higher value !
// #define  DEBUG_FIELD 1    

template
<class Stepper>
class TIntegrationDriver // : public GUIntegrationDriver
{ 
    //  Stepsize can increase by no more than 5.0
    //           and decrease by no more than 1/10. = 0.1
    //
    static const double max_step_increase = 5.0;
    static const double max_step_decrease = 0.1;

    //  The (default) maximum number of steps is Base
    //  divided by the order of Stepper
    //
    static const int  fMaxStepBase = 250;  // Was 5000

    public:  // with description

    typedef Stepper T_Stepper;
    
    static const size_t Nvar = Stepper::Nvar;
    typedef StaticVector<double, Nvar> StateVec;
    typedef BlazePairType<Nvar>        BlazeOutVec; 

    inline bool
        AccurateAdvance(GUFieldTrack& y_current,
                double     hstep,
                double     eps,
                double hinitial=0.0 );

    // ---------------------------------------------------------
    //----------------------------------------------------------------------

    bool  QuickAdvance(       
            GUFieldTrack& y_posvel,         // INOUT
            const StateVec &dydx,  
            double     hstep,       // In
            double&    dchord_step,
            double&    dyerr );

    // --------------------------------------------------------------------------

    //  Constructor
    //
    TIntegrationDriver( double                hminimum, 
            T_Stepper *pStepper,
            int                   numComponents=6,
            int                   statisticsVerbose=1);

    // ---------------------------------------------------------
    //  Destructor
    //
    ~TIntegrationDriver()
    { 
        if( fStatisticsVerboseLevel > 1 )
        {
            PrintStatisticsReport();
        }
    }


    // ---------------------------------------------------------

    inline double GetHmin()   const { return fMinimumStep; } 
    inline double Hmin()      const { return fMinimumStep; }
    inline double GetSafety() const { return safety; }
    inline double GetPowerShrink() const { return fPowerShrink; } 
    inline double GetPowerGrow()  const { return fPowerGrow; }
    inline double GetErrcon() const { return errcon; }
    inline void SetHmin(double newval) { fMinimumStep = newval;         } 

    inline
        double ComputeAndSetErrcon()
        {
            errcon = std::pow(max_step_increase/safety,1.0/fPowerGrow);
            return errcon;
        } 

    inline
        void ReSetParameters(double new_safety=0.9)
        {
            safety = new_safety;
            fPowerShrink = -1.0 /
                    pIntStepper->T_Stepper::IntegratorOrder();
            fPowerGrow  = -1.0 / (1.0 + 
                    pIntStepper->T_Stepper::IntegratorOrder());
            ComputeAndSetErrcon();
        }

    inline void SetSafety(double val) { safety=val;     ComputeAndSetErrcon(); }
    inline void SetPowerGrow(double  val) { fPowerGrow=val;     ComputeAndSetErrcon();  }
    inline void SetErrcon(double val) { errcon=val; }

    inline void RenewStepperAndAdjust(T_Stepper *pItsStepper) 
        {  
            pIntStepper = pItsStepper; 
            ReSetParameters();
        }

    inline const T_Stepper* GetStepper() const { return pIntStepper; }
    inline  T_Stepper* GetStepper()  { return pIntStepper; }
    inline int GetMaxNoSteps() const { return fMaxNoSteps; }
    inline void SetMaxNoSteps(int val) { fMaxNoSteps= val; }

    inline
        StateVec GetDerivatives(const GUFieldTrack &y_curr)
        { 
            double  tmpValArr[GUFieldTrack::ncompSVEC];
            y_curr.DumpToArray(tmpValArr);
            StateVec tmpValArrv(Nvar, tmpValArr);
            return pIntStepper->
                T_Stepper::RightHandSide(tmpValArrv);
        }

    inline double GetVerboseLevel() const {    return fVerboseLevel;  } 
    inline   void SetVerboseLevel(int newLevel) { fVerboseLevel= newLevel; }
    inline double GetSmallestFraction() const { return fSmallestFraction; } 

    void
        OneGoodStep(StateVec &y,        // InOut
                const StateVec &dydx,
                double& x,         // InOut
                double htry,
                double eps_rel_max,
                double& hdid,      // Out
                double& hnext );   // Out

    // ---------------------------------------------------------------------------



    public:  // without description


    void SetSmallestFraction(double newFraction)
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

    protected:  // without description
    void
        WarnSmallStepSize( double hnext, double hstep, 
                double h, double xDone,
                int nstp);

    // ---------------------------------------------------------
    void
        WarnTooManyStep( double x1start, 
                double x2end, 
                double xCurrent)
        {
            std::ostringstream message;
            message << "The number of steps used in the Integration driver"
                << " (Runge-Kutta) is too many." << std::endl
                << "Integration of the interval was not completed !" << std::endl
                << "Only a " << (xCurrent-x1start)*100/(x2end-x1start)
                << " % fraction of it was done.";
            G4Exception("WarnTooManyStep()", "GeomField1001",
                    JustWarning, message);
        }

    // ---------------------------------------------------------
    void
        WarnEndPointTooFar (double endPointDist, 
                double   h , 
                double  eps,
                int     dbg)
        {
            static G4ThreadLocal double maxRelError=0.0;
            bool isNewMax, prNewMax;

            isNewMax = endPointDist > (1.0 + maxRelError) * h;
            prNewMax = endPointDist > (1.0 + 1.05 * maxRelError) * h;
            if( isNewMax ) { maxRelError= endPointDist / h - 1.0; }

            if( dbg && (h > G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) 
                    && ( (dbg>1) || prNewMax || (endPointDist >= h*(1.+eps) ) ) )
            { 
                static G4ThreadLocal int noWarnings = 0;
                std::ostringstream message;
                if( (noWarnings ++ < 10) || (dbg>2) )
                {
                    message << "The integration produced an end-point which " << std::endl
                        << "is further from the start-point than the curve length."
                        << std::endl;
                }
                message << "  Distance of endpoints = " << endPointDist
                    << ", curve length = " << h << std::endl
                    << "  Difference (curveLen-endpDist)= " << (h - endPointDist)
                    << ", relative = " << (h-endPointDist) / h 
                    << ", epsilon =  " << eps;
                G4Exception("WarnEndPointTooFar()", "GeomField1001",
                        JustWarning, message);
            }
        }

    // ---------------------------------------------------------
    // ---------------------------------------------------------------------------

    void PrintStatus( const double*   StartArr,  
            double          xstart,
            const StateVec   &CurrentArr, 
            double          xcurrent,
            double          requestStep, 
            int             subStepNo)
        // Potentially add as arguments:  
        //                                 <dydx>           - as Initial Force
        //                                 stepTaken(hdid)  - last step taken
        //                                 nextStep (hnext) - proposal for size
    {
        GUFieldTrack  StartFT(ThreeVector(0,0,0),
                ThreeVector(0,0,0), 0., 0., 0., 0. );
        GUFieldTrack  CurrentFT (StartFT);

        StartFT.LoadFromArray( StartArr, fNoIntegrationVariables); 
        StartFT.SetCurveLength( xstart);
        double CurrentArrv[GUFieldTrack::ncompSVEC];
        for(size_t i = 0; i < fNoVars; i ++) 
        {CurrentArrv[i] = CurrentArr[i];}
        CurrentFT.LoadFromArray( CurrentArrv, fNoIntegrationVariables); 
        CurrentFT.SetCurveLength( xcurrent );

        PrintStatus(StartFT, CurrentFT, requestStep, subStepNo ); 
    }

    // ---------------------------------------------------------------------------


template
<class Stepper>
    void PrintStatus(
            const GUFieldTrack&  StartFT,
            const GUFieldTrack&  CurrentFT, 
            double             requestStep, 
            int                subStepNo)
    {
        int verboseLevel= fVerboseLevel;
        static G4ThreadLocal int noPrecision= 5;
        int oldPrec= std::cout.precision(noPrecision);
        // std::cout.setf(ios_base::fixed,ios_base::floatfield);

        const ThreeVector StartPosition=       StartFT.GetPosition();
        const ThreeVector StartMomentumDir=   StartFT.GetMomentumDir();
        const ThreeVector CurrentPosition=     CurrentFT.GetPosition();
        const ThreeVector CurrentMomentumDir= CurrentFT.GetMomentumDir();

        double  DotStartCurrentVeloc= StartMomentumDir.dot(CurrentMomentumDir);

        double step_len= CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
        double subStepSize = step_len;

        if( (subStepNo <= 1) || (verboseLevel > 3) )
        {
            subStepNo = - subStepNo;        // To allow printing banner

            std::cout << std::setw( 6)  << " " << std::setw( 25)
                << " TIntegrationDriver: Current Position  and  Direction" << " "
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
            //    << " and MomentumDir is " << CurrentMomentumDir << std::endl;
            // std::cout << "Step taken was " << step_len  
            //    << " out of PhysicalStep= " <<  requestStep << std::endl;
            // std::cout << "Final safety is: " << safety << std::endl;
            // std::cout << "Chord length = " << (CurrentPosition-StartPosition).mag()
            //        << std::endl << std::endl; 
        }
        std::cout.precision(oldPrec);
    }

    // ---------------------------------------------------------------------------
    void PrintStat_Aux(
            const GUFieldTrack&  aFieldTrack,
            double             requestStep, 
            double             step_len,
            int                subStepNo,
            double             subStepSize,
            double             dotVeloc_StartCurr);

    // ---------------------------------------------------------------------------
    void PrintStatisticsReport();


    // ---------------------------------------------------------------------------


#ifdef QUICK_ADV_TWO
    bool  QuickAdvance(       
            double     yarrin[],    // In
            const StateVec     &dydx,  
            double     hstep,       // In
            double     yarrout[],
            double&    dchord_step,
            double&    dyerr )      // In length
    {
        G4Exception("QuickAdvance()", "GeomField0001",
                FatalException, "Not yet implemented.");
        dyerr = dchord_step = hstep * yarrin[0] * dydx[0];
        yarrout[0]= yarrin[0];
    }
#endif

    private:

    TIntegrationDriver(const TIntegrationDriver&);
    TIntegrationDriver& operator=(const TIntegrationDriver&);
    // Private copy constructor and assignment operator.

    private:

    // ---------------------------------------------------------------
    //  INVARIANTS 

    double  fMinimumStep;
    // Minimum Step allowed in a Step (in absolute units)
    double  fSmallestFraction;      //   Expected range 1e-12 to 5e-15;  
    // Smallest fraction of (existing) curve length - in relative units
    //  below this fraction the current step will be the last 

    const int  fNoIntegrationVariables;  // Number of Variables in integration
    const int  fMinNoVars;               // Minimum number for GUFieldTrack
    const int  fNoVars;                  // Full number of variable

    int   fMaxNoSteps;

    double safety;
    double fPowerShrink;   //  exponent for shrinking
    double fPowerGrow;    //  exponent for growth
    double errcon;
    // Parameters used to grow and shrink trial stepsize.

    int    fStatisticsVerboseLevel;

    // ---------------------------------------------------------------
    // DEPENDENT Objects
    //T_Stepper 
    T_Stepper *pIntStepper;

    // ---------------------------------------------------------------
    //  STATE

    int  fNoTotalSteps, fNoBadSteps, fNoSmallSteps, fNoInitialSmallSteps; 
    double fDyerr_max, fDyerr_mx2;
    double fDyerrPos_smTot, fDyerrPos_lgTot, fDyerrVel_lgTot; 
    double fSumH_sm, fSumH_lg; 
    // Step Statistics 

    int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
    // Could be varied during tracking - to help identify issues

};

class TIntegrationDriver::
  TIntegrationDriver( double hminimum, 
                Stepper *pStepper,
                int                   numComponents=6,
                int                   statisticsVerbose=1)
        : fSmallestFraction( 1.0e-12 ), 
        fNoIntegrationVariables(numComponents), 
        fMinNoVars(12), 
        fNoVars( std::max( fNoIntegrationVariables, fMinNoVars )),
        fStatisticsVerboseLevel(statisticsVerbose),
        fNoTotalSteps(0),  fNoBadSteps(0), fNoSmallSteps(0),
        fNoInitialSmallSteps(0), 
        fDyerr_max(0.0), fDyerr_mx2(0.0), 
        fDyerrPos_smTot(0.0), fDyerrPos_lgTot(0.0), fDyerrVel_lgTot(0.0), 
        fSumH_sm(0.0), fSumH_lg(0.0),
        fVerboseLevel(0), 
        GUIntegrationDriver(hminimum, pStepper, numComponents, statisticsVerbose)
{  
    // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
    // is required. For proper time of flight and spin,  fMinNoVars must be 12

    RenewStepperAndAdjust( pStepper );
    MinimumStep= hminimum;
    fMaxNoSteps = fMaxStepBase / 
    pIntStepper->T_Stepper::IntegratorOrder();
#ifdef DEBUG_FIELD
    fVerboseLevel=2;
#endif

    if( (fVerboseLevel > 0) || (fStatisticsVerboseLevel > 1) )
    {
            std::cout << "MagIntDriver version: Accur-Adv: "
                << "invE_nS, QuickAdv-2sqrt with Statistics "
#ifdef FIELD_STATS
                << " enabled "
#else
                << " disabled "
#endif
                << std::endl;
    }
}

void
OneGoodStep(StateVec &y,        // InOut
        const StateVec &dydx,
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

    BlazeOutVec yVec;
    h = htry ; // Set stepsize to the initial trial value

    double inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

    double errpos_sq=0.0;    // square of displacement error
    double errvel_sq=0.0;    // square of momentum vector difference
    double errspin_sq=0.0;   // square of spin vector difference

    int iter;

    static G4ThreadLocal int tot_no_trials=0; 
    const int max_trials=100; 

    ThreeVector Spin(y[9],y[10],y[11]);
    bool     hasSpin= (Spin.mag2() > 0.0); 

    for (iter=0; iter<max_trials ;iter++)
    {
        tot_no_trials++;
        yVec = pIntStepper->T_Stepper::StepWithError(y,dydx,h); 
        //            *******

        double eps_pos = eps_rel_max * std::max(h, fMinimumStep); 
        double inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos); 

        // Evaluate accuracy
        //
        errpos_sq =  sqr(yVec.err[0]) + sqr(yVec.err[1]) + sqr(yVec.err[2]) ;
        errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

        // Accuracy for momentum
        errvel_sq =  (sqr(yVec.err[3]) + sqr(yVec.err[4]) + sqr(yVec.err[5]) )
            / (sqr(y[3]) + sqr(y[4]) + sqr(y[5]) );
        errvel_sq *= inv_eps_vel_sq;
        errmax_sq = std::max( errpos_sq, errvel_sq ); // Square of maximum error

        if( hasSpin )
        { 
            // Accuracy for spin
            errspin_sq =  ( sqr(yVec.err[9]) + sqr(yVec.err[10]) + sqr(yVec.err[11]) )
                /  ( sqr(y[9]) + sqr(y[10]) + sqr(y[11]) );
            errspin_sq *= inv_eps_vel_sq;
            errmax_sq = std::max( errmax_sq, errspin_sq ); 
        }

        if ( errmax_sq <= 1.0 )  { break; } // Step succeeded. 

        // Step failed; compute the size of retrial Step.
        htemp = GetSafety()*h* G4Pow::GetInstance()->G4Pow::powA( errmax_sq, 0.5*GetPowerShrink() );

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

#ifdef FIELD_STATS
    // Sum of squares of position error // and momentum dir (underestimated)
    fSumH_lg += h; 
    fDyerrPos_lgTot += errpos_sq;
    fDyerrVel_lgTot += errvel_sq * h * h; 
#endif

    // Compute size of next Step
    if (errmax_sq > errcon*errcon)
    { 
        hnext = GetSafety()*h*G4Pow::GetInstance()->G4Pow::powA(errmax_sq, 0.5*GetPowerGrow());
    }
    else
    {
        hnext = max_step_increase*h ; // No more than a factor of 5 increase
    }
    x += (hdid = h);

    y = yVec.out;

    return;
}   // end of  OneGoodStep .............................

//----------------------------------------------------------------------

//  This method computes new step sizes - but does not limit changes to
//   within  certain factors
// 
double 
template <class Stepper>
TIntegrationDriver::ComputeNewStepSize( 
        double  errMaxNorm,    // max error  (normalised)
        double  hstepCurrent)  // current step size
{
    double hnew;

    // Compute size of next Step for a failed step
    if(errMaxNorm > 1.0 )
    {
        // Step failed; compute the size of retrial Step.
        hnew = GetSafety()*hstepCurrent*G4Pow::GetInstance()->G4Pow::powA(errMaxNorm,GetPowerShrink()) ;
    } else if(errMaxNorm > 0.0 ) {
        // Compute size of next Step for a successful step
        hnew = GetSafety()*hstepCurrent*G4Pow::GetInstance()->G4Pow::powA(errMaxNorm,GetPowerGrow()) ;
    } else {
        // if error estimate is zero (possible) or negative (dubious)
        hnew = max_step_increase * hstepCurrent; 
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
template <class Stepper>
TIntegrationDriver::ComputeNewStepSize_WithinLimits( 
        double  errMaxNorm,    // max error  (normalised)
        double  hstepCurrent)  // current step size
{
    double hnew;

    // Compute size of next Step for a failed step
    if (errMaxNorm > 1.0 )
    {
        // Step failed; compute the size of retrial Step.
        hnew = GetSafety()*hstepCurrent*G4Pow::GetInstance()->G4Pow::powA(errMaxNorm,GetPowerShrink()) ;

        if (hnew < max_step_decrease*hstepCurrent)
        {
            hnew = max_step_decrease*hstepCurrent ;
            // reduce stepsize, but no more
            // than this factor (value= 1/10)
        }
    }
    else
    {
        // Compute size of next Step for a successful step
        if (errMaxNorm > errcon)
        { hnew = GetSafety()*hstepCurrent*G4Pow::GetInstance()->G4Pow::powA(errMaxNorm,GetPowerGrow()); }
        else  // No more than a factor of 5 increase
        { hnew = max_step_increase * hstepCurrent; }
    }
    return hnew;
}

bool
template
<class Stepper>
class TIntegrationDriver::AccurateAdvance(GUFieldTrack& y_current,
        double     hstep,
        double     eps,
        double hinitial=0.0 )
{
    // Runge-Kutta driver with adaptive stepsize control. Integrate starting
    // values at y_current over hstep x2 with accuracy eps. 
    // On output ystart is replaced by values at the end of the integration 
    // interval. RightHandSide is the right-hand side of ODE system. 
    // The source is similar to odeint routine from NRC p.721-722 .

    int nstp, i, no_warnings=0;
    double x, hnext, hdid, h;

#ifdef DEBUG_FIELD
    static int dbg=1;
    static int nStpPr=50;   // For debug printing of long integrations
    double ySubStepStart[GUFieldTrack::ncompSVEC];
    GUFieldTrack  yFldTrkStart(y_current);
#endif

    StateVec y; 
    StateVec dydx;
    double ystart[GUFieldTrack::ncompSVEC], yEnd[GUFieldTrack::ncompSVEC]; 
    double  x1, x2;
    bool succeeded = true, lastStepSucceeded;

    double startCurveLength;

    int  noFullIntegr=0, noSmallIntegr = 0 ;
    static G4ThreadLocal int  noGoodSteps =0 ;  // Bad = chord > curve-len 

    GUFieldTrack yStartFT(y_current);

    //  Ensure that hstep > 0
    //
    if( hstep <= 0.0 )
    { 
        if(hstep==0.0)
        {
            std::ostringstream message;
            message << "Proposed step is zero; hstep = " << hstep << " !";
            G4Exception("AccurateAdvance()", 
                    "GeomField1001", JustWarning, message);
            return succeeded; 
        }
        else
        { 
            std::ostringstream message;
            message << "Invalid run condition." << std::endl
                << "Proposed step is negative; hstep = " << hstep << "." << std::endl
                << "Requested step cannot be negative! Aborting event.";
            G4Exception("AccurateAdvance()", 
                    "GeomField0003", EventMustBeAborted, message);
            return false;
        }
    }

    y_current.DumpToArray( ystart );

    startCurveLength= y_current.GetCurveLength();
    x1= startCurveLength; 
    x2= x1 + hstep;

    if ( (hinitial > 0.0) && (hinitial < hstep)
            && (hinitial > perMillion * hstep) ){
        h = hinitial;
    }else{  //  Initial Step size "h" defaults to the full interval
        h = hstep;
    }

    x = x1;

    for (i=0;i<fNoVars;i++)  { y[i] = ystart[i]; }

    bool lastStep= false;
    nstp=1;

    do
    {
        ThreeVector StartPos( y[0], y[1], y[2] );

#ifdef DEBUG_FIELD
        double xSubStepStart= x; 
        for (i=0;i<fNoVars;i++)  { ySubStepStart[i] = y[i]; }
        yFldTrkStart.LoadFromArray(ySubStepStart, fNoIntegrationVariables);
        yFldTrkStart.SetCurveLength(x);
#endif

        // Old method - inline call to Equation of Motion
        //   pIntStepper->T_Stepper::RightHandSide( y, dydx );
        // New method allows to cache field, or state (eg momentum magnitude)
        dydx = pIntStepper->
            T_Stepper::RightHandSide( y );

        fNoTotalSteps++;

        // Perform the Integration
        //      
        if( h > fMinimumStep )
        { 
            OneGoodStep(y,dydx,x,h,eps,hdid,hnext) ;
            //--------------------------------------
            lastStepSucceeded= (hdid == h);   
#ifdef DEBUG_FIELD
            if (dbg>2) {
                PrintStatus( ySubStepStart, xSubStepStart, y, x, h,  nstp); // Only
            }
#endif
        }
        else
        {
            GUFieldTrack yFldTrk( ThreeVector(0,0,0), 
                    ThreeVector(0,0,0), 0., 0., 0., 0. );
            double dchord_step, dyerr, dyerr_len;   // What to do with these ?
            double yv[GUFieldTrack::ncompSVEC];
            for(size_t i = 0; i < Nvar; i ++){ yv[i] = y[i];}
            yFldTrk.LoadFromArray(yv, fNoIntegrationVariables); 
            yFldTrk.SetCurveLength( x );

            QuickAdvance( yFldTrk, dydx, h, dchord_step, dyerr_len ); 
            //-----------------------------------------------------
            yFldTrk.DumpToArray(yv);
            for(size_t i = 0; i < Nvar; i ++){ y[i] = yv[i];}

#ifdef FIELD_STATS
            fNoSmallSteps++; 
            if ( dyerr_len > fDyerr_max)  { fDyerr_max= dyerr_len; }
            fDyerrPos_smTot += dyerr_len;
            fSumH_sm += h;  // Length total for 'small' steps
            if (nstp<=1)  { fNoInitialSmallSteps++; }
#endif
#ifdef DEBUG_FIELD
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
            if( h == 0.0 )
            { 
                G4Exception("AccurateAdvance()",
                        "GeomField0003", FatalException,
                        "Integration Step became Zero!"); 
            }
            dyerr = dyerr_len / h;
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

#ifdef  DEBUG_FIELD
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
        double endPointDist= (EndPos-StartPos).mag(); 
        if ( endPointDist >= hdid*(1.+perMillion) )
        {
            fNoBadSteps++;

            // Issue a warning only for gross differences -
            // we understand how small difference occur.
            if ( endPointDist >= hdid*(1.+perThousand) )
            { 
#ifdef DEBUG_FIELD
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
            if(std::fabs(hnext) <= Hmin())
            {
#ifdef  DEBUG_FIELD
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
                h = Hmin();
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
#ifdef DEBUG_FIELD
                if (dbg>2)
                {
                    int prec= std::cout.precision(12); 
                    std::cout << "Warning: TMagIntegratorDriver::AccurateAdvance"
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

    for (i=0;i<fNoVars;i++)  { yEnd[i] = y[i]; }

    // Put back the values.
    y_current.LoadFromArray( yEnd, fNoIntegrationVariables );
    y_current.SetCurveLength( x );

    if(nstp > fMaxNoSteps)
    {
        no_warnings++;
        succeeded = false;
#ifdef DEBUG_FIELD
        if (dbg)
        {
            WarnTooManyStep( x1, x2, x );  //  Issue WARNING
            PrintStatus( yEnd, x1, y, x, hstep, -nstp);
        }
#endif
    }

#ifdef DEBUG_FIELD
    if( dbg && no_warnings )
    {
        std::cerr << "TMagIntegratorDriver exit status: no-steps " << nstp <<std::endl;
        PrintStatus( yEnd, x1, y, x, hstep, nstp);
    }
#endif

    return succeeded;
}  // end of AccurateAdvance ...........................


bool
template <class Stepper> TIntegrationDriver::
QuickAdvance(       
            GUFieldTrack& y_posvel,         // INOUT
            const StateVec &dydx,  
            double     hstep,       // In
            double&    dchord_step,
            double&    dyerr )
{
    double dyerr_pos_sq, dyerr_mom_rel_sq;  
    double s_start;
    double dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

    static G4ThreadLocal int no_call=0; 
    no_call ++; 

    // Move data into array
    double yarrin[GUFieldTrack::ncompSVEC]; 
    y_posvel.DumpToArray( yarrin );      //  yarrin  <== y_posvel 
    s_start = y_posvel.GetCurveLength();

    StateVec yarrinv(Nvar, yarrin);
    // Do an Integration Step
    BlazeOutVec yVec(pIntStepper->
                     T_Stepper::StepWithError(yarrinv, dydx, hstep)); 
    //            *******

    // Estimate curve-chord distance
    dchord_step= pIntStepper->T_Stepper:: DistChord();
    //                         *********

    // Put back the values.  yarrout ==> y_posvel
    double yarrout[GUFieldTrack::ncompSVEC];
    for(size_t i = 0; i < Nvar; i ++){ yarrout[i] = yVec.out[i]; }
        y_posvel.LoadFromArray( yarrout, fNoIntegrationVariables );
    y_posvel.SetCurveLength( s_start + hstep );

#ifdef  DEBUG_FIELD
    if(fVerboseLevel>2)
    {
        std::cout << "G4MagIntDrv: Quick Advance" << std::endl;
        PrintStatus( yarrin, s_start, yVec.out, s_start+hstep, hstep,  1); 
    }
#endif

    // A single measure of the error   
    //      TO-DO :  account for  energy,  spin, ... ? 
    vel_mag_sq   = ( sqr(yVec.out[3])+sqr(yVec.out[4])+sqr(yVec.out[5]) );
    inv_vel_mag_sq = 1.0 / vel_mag_sq; 
    dyerr_pos_sq = ( sqr(yVec.err[0])+sqr(yVec.err[1])+sqr(yVec.err[2]));
    dyerr_mom_sq = ( sqr(yVec.err[3])+sqr(yVec.err[4])+sqr(yVec.err[5]));
    dyerr_mom_rel_sq = dyerr_mom_sq * inv_vel_mag_sq;

    // Calculate also the change in the momentum squared also ???
    // double veloc_square = y_posvel.GetVelocity().mag2();
    // ...

#ifdef RETURN_A_NEW_STEP_LENGTH
    // The following step cannot be done here because "eps" is not known.
    dyerr_len = 1.0/vdt::fast_isqrt_general( dyerr_len_sq, 2); 
    dyerr_len_sq /= eps ;

    // Look at the velocity deviation ?
    //  sqr(yerr_vec[3])+sqr(yerr_vec[4])+sqr(yerr_vec[5]));

    // Set suggested new step
    hstep= ComputeNewStepSize( dyerr_len, hstep);
#endif

    if( dyerr_pos_sq > ( dyerr_mom_rel_sq * sqr(hstep) ) )
    {
        dyerr = 1.0/vdt::fast_isqrt_general(dyerr_pos_sq, 2);
    }
    else
    {
        // Scale it to the current step size - for now
        dyerr = 1.0/vdt::fast_isqrt_general(dyerr_mom_rel_sq, 2) * hstep;
    }

    return true;
}

void 
template <class Stepper> TIntegrationDriver::
PrintStatisticsReport()
{
    int noPrecBig= 6;
    int oldPrec= std::cout.precision(noPrecBig);

    std::cout << "TIntegrationDriver Statistics of steps undertaken. " << std::endl;
    std::cout << "TIntegrationDriver: Number of Steps: "
        << " Total= " <<  fNoTotalSteps
        << " Bad= "   <<  fNoBadSteps 
        << " Small= " <<  fNoSmallSteps 
        << " Non-initial small= " << (fNoSmallSteps-fNoInitialSmallSteps)
        << std::endl;

#ifdef FIELD_STATS
    std::cout << "MID dyerr: " 
        << " maximum= " << fDyerr_max 
        << " Sum small= " << fDyerrPos_smTot 
        << " 1.0/vdt::fast_isqrt_general(Sum large^2): pos= " << 1.0/vdt::fast_isqrt_general(fDyerrPos_lgTot, 2)
        << " vel= " << 1.0/vdt::fast_isqrt_general( fDyerrVel_lgTot, 2)
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

void
template <class Stepper> TIntegrationDriver::
        WarnSmallStepSize( double hnext, double hstep, 
                double h, double xDone,
                int nstp)
        {
            static G4ThreadLocal int noWarningsIssued =0;
            const  int maxNoWarnings =  10;   // Number of verbose warnings
            std::ostringstream message;
            if( (noWarningsIssued < maxNoWarnings) || fVerboseLevel > 10 )
            {
                message << "The stepsize for the next iteration, " << hnext
                    << ", is too small - in Step number " << nstp << "." << std::endl
                    << "The minimum for the driver is " << Hmin()  << std::endl
                    << "Requested integr. length was " << hstep << " ." << std::endl
                    << "The size of this sub-step was " << h     << " ." << std::endl
                    << "The integrations has already gone " << xDone;
            }
            else
            {
                message << "Too small 'next' step " << hnext
                    << ", step-no: " << nstp << std::endl
                    << ", this sub-step: " << h     
                    << ",  req_tot_len: " << hstep 
                    << ", done: " << xDone << ", min: " << Hmin();
            }
            G4Exception("WarnSmallStepSize()", "GeomField1001",
                    JustWarning, message);
            noWarningsIssued++;
        }

void
template <class Stepper> TIntegrationDriver::
 PrintStat_Aux(
        const GUFieldTrack&  aFieldTrack,
        double             requestStep, 
        double             step_len,
        int                subStepNo,
        double             subStepSize,
        double             dotVeloc_StartCurr)
{
    const ThreeVector Position=      aFieldTrack.GetPosition();
    const ThreeVector MomentumDir=   aFieldTrack.GetMomentumDir();

    if( subStepNo >= 0)
    {
        std::cout << std::setw( 5) << subStepNo << " ";
    }
    else
    {
        std::cout << std::setw( 5) << "Start" << " ";
    }
    double curveLen= aFieldTrack.GetCurveLength();
    std::cout << std::setw( 7) << curveLen;
    std::cout << std::setw( 9) << Position.x() << " "
        << std::setw( 9) << Position.y() << " "
        << std::setw( 9) << Position.z() << " "
        << std::setw( 8) << MomentumDir.x() << " "
        << std::setw( 8) << MomentumDir.y() << " "
        << std::setw( 8) << MomentumDir.z() << " ";
    int oldprec= std::cout.precision(3);
    std::cout << std::setw( 8) << MomentumDir.mag2()-1.0 << " ";
    std::cout.precision(6);
    std::cout << std::setw(10) << dotVeloc_StartCurr << " ";
    std::cout.precision(oldprec);
    std::cout << std::setw( 7) << aFieldTrack.GetKineticEnergy();
    std::cout << std::setw(12) << step_len << " ";

    static G4ThreadLocal double oldCurveLength= 0.0;
    static G4ThreadLocal double oldSubStepLength= 0.0;
    static G4ThreadLocal int oldSubStepNo= -1;

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


#endif /* TIntegrationDriver_Def */
