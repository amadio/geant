//
//
// Class description:
//
//  'Drives' integration in multiple lanes, asking the Integrator Stepper to do
//     the integration work for a few lanes at a time.   Similarly to the original
//     sequential version, it ensures that the step's error remains within acceptable
//     bounds.
//   Different approaches are available to handle 'divergence' - i.e. some lanes 
//     finished while others are continuing (failed attempts or require further
//     sub-steps):
//     - OneStep() just finishes the integration of all vector lanes (simplest)
//     - KeepStepping() stores the results for finished lanes, and reuses its
//         working array(s) for further lanes  (optimisation)
//
//   Started on GUVIntegrationDriver; introduces template on backend.
//   First implementation:  Ananya.   Supervision:  J. Apostolakis
//
// History:
// - First version created by Ananya,  Feb 2016, 
//   Revisions by J.Apostolakis,   Apr-Jun 2016
// -----------------------------------------------------------------------------

#ifndef TemplateGUIntegrationDriver_Def
#define TemplateGUIntegrationDriver_Def

#include "TemplateGUFieldTrack.h"
#include "AlignedBase.h"
#include "FieldTrack.h"

#include "TemplateGUVIntegrationStepper.h"

#include "base/Vector.h"

// Needed for scalar stepper for new constructor (KeepStepping)
#include "GUVIntegrationStepper.h"

// Needed to send work to scalar driver to deal with 1/2 remaining lanes
#include "GUIntegrationDriver.h"
#include "GUFieldTrack.h"

#define NEWACCURATEADVANCE

template <class Backend>
class TemplateGUIntegrationDriver : public AlignedBase
{
   public:  // with description

     typedef                   typename Backend::precision_v  Double_v;
     typedef                   typename Backend::bool_v       Bool_v;
     typedef vecgeom::Vector3D<typename Backend::precision_v> ThreeVectorSimd; 

     TemplateGUIntegrationDriver( double hminimum,  //same 
                                  TemplateGUVIntegrationStepper<Backend> *pStepper,
                                  int    numberOfComponents  = 6,
                                  int    statisticsVerbosity = 1                  );

     TemplateGUIntegrationDriver( const TemplateGUIntegrationDriver& );
       // Copy constructor used to create Clone method
     ~TemplateGUIntegrationDriver();

     Bool_v  AccurateAdvance( const TemplateGUFieldTrack<Backend>& y_current,
                              const Double_v  charge,
                                    Double_v  hstep,
                                    double    eps,      //same      // Requested y_err/hstep
                                    TemplateGUFieldTrack<Backend>& yOutput   );                            
                                    // Double_v  hinitial=0.0                  );  // Suggested 1st interval
       // Above drivers for integrator (Runge-Kutta) with stepsize control. 
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps. 
       // On output ystart is replaced by value at end of interval. 
       // The concept is similar to the odeint routine from NRC p.721     
#ifdef NEWACCURATEADVANCE
     // 2nd AccurateAdvance
     // invovles track insertion etc 
     // succeeded[] of length nTracks
     void  AccurateAdvance( /*const*/ FieldTrack yInput[],
                            const double      charge[],
                                  double      hstep[],
                                  double      epsilon,
                                  FieldTrack  yOutput[],
                                  int         nTracks, 
                                  bool        succeeded[] );

     void  OneStep(       Double_v  ystart[], // Like old RKF45step()
                    const Double_v  dydx[],
                    const Double_v  charge,
                    const Double_v  epsRel,   // constant ? -- ok to have it per track! 
                          Double_v& x,
                          Double_v  htry,
                          Double_v& hdid,
                          Double_v& hnext ) ;

     void  InitializeLanes( /*const*/ FieldTrack  yInput[],
                                const double      hstep [],
                                const double      chargeIn[],
                                      Double_v    yLane[],
                                      Double_v &  hStepLane,
                                      Double_v &  chargeLane,                                                 
                                      Double_v &  startCurveLength);

     // Returns isDoneThisLane value i.e. whether the lane is done or not
     // True if the last lane left is an h<=0 lane, false for all else
     bool  InsertNewTrack( /*const*/ FieldTrack yInput[],
                            const double     hstep[],
                            const double     chargeIn[],
                            const int        currIndex,
                                  int&       trackNextInput,
                                  bool       succeeded[],
                                  Double_v   y[],
                                  Double_v  &hStepLane,
                                  Double_v  &chargeLane,                                               
                                  Double_v  &startCurveLength );

     void  StoreOutput( const Double_v   y[],
                        const Double_v   x,
                              FieldTrack yOutput[],
                              int        currIndex,
                              double     hstep[],
                              bool       succeeded[]);

     void  SetNTracks( int nTracks );

     void  KeepStepping(       Double_v  ystart[], // Like old RKF45step()
                               Double_v  dydx[],
                               Double_v& x,                               
                         const Double_v  charge,
                         const Double_v  htry,
                         const Double_v  eps,      //  memb variables ?
                         const Double_v  hTotalLane,
                      //       Double_v& hdid,
                               Double_v& hnext,
                               Double_v& hTotalDoneSoFar,
                               Double_v& stepsLane    // Number of steps 
        ) ;

     TemplateGUIntegrationDriver( double hminimum,  //same 
                                  TemplateGUVIntegrationStepper<Backend> *pStepper,
                                  GUVIntegrationStepper                  *pScalarStepper,
                                  int    numberOfComponents  = 6,
                                  int    statisticsVerbosity = 1                        );

     TemplateGUIntegrationDriver( double hminimum,  //same 
                                  TemplateGUVIntegrationStepper<Backend> *pStepper,
                                  GUVIntegrationStepper                  *pScalarStepper,
                                  GUIntegrationDriver                    *pScalarDriver,
                                  int    numberOfComponents  = 6,
                                  int    statisticsVerbosity = 1                        );

     void SetPartDebug (bool debugValue);

     // true for OneStep, false for KeepStepping
     void SetSteppingMethod (bool steppingMethod);
#endif

     Bool_v  QuickAdvance( TemplateGUFieldTrack<Backend>& y_posvel,        // INOUT
                           const Double_v      dydx[],
                           const Double_v      charge,          // In 
                                 Double_v      hstep,           // IN
      #ifdef USE_DCHORD
                                 Double_v&     dchord_step, //take out
      #endif                       
                                 Double_v&     dyerr_pos_sq,
                                 Double_v&     dyerr_mom_rel_sq ) ;
       // New QuickAdvance that also just tries one Step
       //    (so also does not ensure accuracy)
       //    but does return the errors in  position and
       //        momentum (normalised: Delta_Integration(p^2)/(p^2) )

     // void  InitializeCharge(Double_v charge) { fpStepper->InitializeCharge(charge);} //remove 
       // Pass needed information and initialize 
     void  DoneIntegration() { fpStepper->GetEquationOfMotion()->InformDone(); } 
       // Pass along information about end of integration - can clears parameters, flag finished

     void ReportStepUnderflow( const std::string method, double x, double xnew, double hOriginal,
                               double hLast, double hNext,
                               double errPosSq, double errMomSq, const double dydx[], int components);

     TemplateGUIntegrationDriver* Clone() const;
       // Create an independent copy of the current object -- including independent 'owned' objects
       // 
       // Question:  If the current object and all sub-objects are const, can it return 'this' ?
     
     TemplateGUVEquationOfMotion<Backend>* GetEquationOfMotion() { return fpStepper->GetEquationOfMotion(); }
     const TemplateGUVEquationOfMotion<Backend>* GetEquationOfMotion() const { return fpStepper->GetEquationOfMotion(); } 

     static constexpr int ncompSVEC = TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC;
     
     // Auxiliary methods
     inline double GetHmin()        const { return fMinimumStep; } 
     inline double GetSafety()      const { return fSafetyFactor; }
     inline double GetPowerShrink() const { return fPowerShrink; }
     inline double GetPowerGrow()   const { return fPowerGrow; } 
     inline double GetErrcon()      const { return fErrcon; }
     
     inline void   GetDerivatives( const TemplateGUFieldTrack<Backend> &y_curr,     // const, INput
                                        Double_v    charge, 
                                        Double_v    dydx[]   );  //       OUTput
        // Accessors.

     inline void RenewStepperAndAdjust(TemplateGUVIntegrationStepper<Backend> *Stepper);
        // Sets a new stepper 'Stepper' for this driver. Then it calls
        // ReSetParameters to reset its parameters accordingly.

     inline void ReSetParameters(double new_safety= 0.9 ); //same 
        //  i) sets the exponents (fPowerGrow & fPowerShrink), 
        //     using the current Stepper's order, 
        // ii) sets the safety
        // ii) calculates "fErrcon" according to the above values.

     inline void SetSafety(double valS) ;       // {fSafetyFactor= valS; }
     inline void SetPowerShrink(double valPs) ; // { fPowerShrink= valPs; } 
     inline void SetPowerGrow (double valPg) ;  // { fPowerGrow= valPg; }
     inline void SetErrcon(double valEc) ;      // { fErrcon= valEc; }
        // When setting safety or fPowerGrow, errcon will be set to a 
        // compatible value.

     inline double ComputeAndSetErrcon();

     inline const TemplateGUVIntegrationStepper<Backend>* GetStepper() const;
     inline TemplateGUVIntegrationStepper<Backend>* GetStepper();

     void  OneGoodStep(       Double_v  ystart[], // Like old RKF45step()
                        const Double_v  dydx[],
                              Double_v  charge,
                              Double_v& x,
                              Double_v  htry,
                              Double_v  eps,      //  memb variables ?
                              Double_v& hdid,
                              Double_v& hnext ) ;
        // This takes one Step that is as large as possible while 
        // satisfying the accuracy criterion of:
        // yerr < eps * |y_end-y_start|

     Double_v ComputeNewStepSize( Double_v  errMaxNorm,    // normalised error
                                  Double_v  hstepCurrent); // current step size
        // Taking the last step's normalised error, calculate
        // a step size for the next step.
        // Do not limit the next step's size within a factor of the
        // current one.

     Double_v ComputeNewStepSize_WithinLimits( Double_v  errMaxNorm,    // normalised error
                                               Double_v  hstepCurrent); // current step size
        // Taking the last step's normalised error, calculate
        // a step size for the next step.
        // Limit the next step's size within a range around the current one.

     inline int    GetMaxNoSteps() const;
     inline void   SetMaxNoSteps( int val); 
        //  Modify and Get the Maximum number of Steps that can be
        //   taken for the integration of a single segment -
        //   (ie a single call to AccurateAdvance).

   public:  // without description

     inline void SetHmin(double newMin) { fMinimumStep= newMin; }
     inline void SetVerboseLevel(int lev) { fVerboseLevel= lev; }
     inline int  GetVerboseLevel() const { return fVerboseLevel; }

     inline double GetSmallestFraction() const { return fSmallestFraction; } 
     void     SetSmallestFraction( double val ); 

   protected:  // without description
     void WarnSmallStepSize( Double_v , // hnext, 
                             Double_v , // hstep, 
                             Double_v , // h,     
                             Double_v , // xDone,
                             int      /*noSteps*/ ) {} //;  //warnings per track, probably neeed to change all to double :ananya
     void WarnTooManySteps ( Double_v , // x1start, 
                             Double_v , // x2end, 
                             Double_v /*xCurrent*/ ) {} //;
     void WarnEndPointTooFar( Double_v  , // endPointDist, 
                              Double_v  , // hStepSize , 
                              Double_v  , // epsilonRelative,
                              int     /*debugFlag*/ ) {} //;
     //  Issue warnings for undesirable situations
     //add index in order to print one at a time :ananya 
     void PrintStatus(  const Double_v*      , // StartArr,
                              Double_v       , // xstart,
                        const Double_v*      , // CurrentArr, 
                              Double_v       , // xcurrent, 
                              Double_v       , // requestStep, 
                              int          /* subStepNo */  ) {} //;
     void PrintStatus(  const TemplateGUFieldTrack<Backend>&  /*StartFT*/,
                        const TemplateGUFieldTrack<Backend>&  /*CurrentFT*/, 
                              double       /*requestStep*/, 
                              int          /*subStepNo*/ ) {} //;
     void PrintStat_Aux( const TemplateGUFieldTrack<Backend>& , // aGUFieldTrack,
                               double      , // requestStep, 
                               double      , // actualStep,
                               int         , // subStepNo,
                               double      , // subStepSize,
                               double      /*dotVelocities*/ ) {} //;       
     //  Verbose output for debugging

     void PrintStatisticsReport() ; // {} //;
     //  Report on the number of steps, maximum errors etc.

#ifdef QUICK_ADV_ARRAY_IN_AND_OUT      
     Bool_v QuickAdvance(      Double_v     yarrin[],     // In
                         const Double_v     dydx[],
                         const Double_v     charge,                               
                               Double_v     hstep,
                               Double_v     yarrout[],    // Out
                               Double_v&    dchord_step,  // Out
                               Double_v&    dyerr );      // in length
#endif

   private:

     TemplateGUIntegrationDriver& operator=(const TemplateGUIntegrationDriver&);
     // Private copy constructor and assignment operator.

   private:

     // ---------------------------------------------------------------
     // DEPENDENT Objects
     TemplateGUVIntegrationStepper<Backend> *fpStepper;
     
     // ---------------------------------------------------------------
     //  INVARIANTS 

     double  fMinimumStep; //same 
     // Minimum Step allowed in a Step (in absolute units)
     double  fSmallestFraction; //same      //   Expected range 1e-12 to 5e-15;  
     // Smallest fraction of (existing) curve length - in relative units
     //  below this fraction the current step will be the last 

     const int  fNoIntegrationVariables;  // Number of Variables in integration
     const int  fMinNoVars;               // Minimum number for TemplateGUFieldTrack<Backend>
     const int  fNoVars;                  // Full number of variable

     int   fMaxNoSteps;
     static const int  fMaxStepBase;  

     double fSafetyFactor;
     double fPowerShrink;   //  exponent for shrinking
     double fPowerGrow;    //  exponent for growth
     double fErrcon;
     // Parameters used to grow and shrink trial stepsize.

     double fSurfaceTolerance; 

     //  Stepsize can increase by no more than 5.0
     //           and decrease by no more than x10. = 0.1 
     static constexpr double fMaxSteppingIncrease = 5.0;
     static constexpr double fMaxSteppingDecrease = 0.1;
     // Maximum stepsize increase/decrease factors.

     int    fStatisticsVerboseLevel;


     // ---------------------------------------------------------------
     //  STATE
    public:
     int    fNoTotalSteps, fNoBadSteps, fNoSmallSteps, fNoInitialSmallSteps; 
     Double_v fDyerr_max, fDyerr_mx2;
     Double_v fDyerrPos_smTot, fDyerrPos_lgTot, fDyerrVel_lgTot; 
     Double_v fSumH_sm, fSumH_lg; 

     int     fOneGoodStepCallsVector, fOneStepCallsScalar;
     // Step Statistics 

     int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
     // Could be varied during tracking - to help identify issues

#ifdef NEWACCURATEADVANCE
     //Variables required for track insertion algorithm
     static constexpr int kVectorSize = vecgeom::kVectorSize; // 4; //can be templated on the backend somehow
     int *fIndex; // or int fIndex[kVectorSize]
     int fNTracks;
     int fStepperCalls = 0;
     GUVIntegrationStepper *fpScalarStepper;
     GUIntegrationDriver *fpScalarDriver;
     bool partDebug = false ; 
     bool oneStep = true; // false for KeepStepping
#endif 
};

// #include "GUIntegratorDriver.icc"

template<class Backend>
constexpr double TemplateGUIntegrationDriver<Backend>::fMaxSteppingIncrease;

template<class Backend>
constexpr double TemplateGUIntegrationDriver<Backend>::fMaxSteppingDecrease;


template <class Backend>
inline
double TemplateGUIntegrationDriver<Backend>
  ::ComputeAndSetErrcon()
{
  fErrcon = std::pow(fMaxSteppingIncrease/fSafetyFactor,1.0/fPowerGrow);
  return fErrcon;
} 

template <class Backend>
inline
void TemplateGUIntegrationDriver<Backend>
  ::ReSetParameters(double new_safety)
{
  fSafetyFactor = new_safety;
  fPowerShrink  = -1.0 / fpStepper->IntegratorOrder();
  fPowerGrow    = -1.0 / (1.0 + fpStepper->IntegratorOrder());
  ComputeAndSetErrcon();
}

template <class Backend>
inline
void TemplateGUIntegrationDriver<Backend>
  ::SetSafety(double val)
{ 
  fSafetyFactor=val;
  ComputeAndSetErrcon();
}

template <class Backend>
inline
void TemplateGUIntegrationDriver<Backend>
  ::SetPowerGrow(double  val)
{ 
  fPowerGrow=val;
  ComputeAndSetErrcon(); 
}

template <class Backend>
inline
void TemplateGUIntegrationDriver<Backend>
  ::SetErrcon(double val)
{ 
  fErrcon=val;
}

template <class Backend>
inline
void TemplateGUIntegrationDriver<Backend>
  ::RenewStepperAndAdjust(TemplateGUVIntegrationStepper<Backend> *pStepper)
{  
  fpStepper = pStepper; 
  ReSetParameters();
}

template <class Backend>
inline
const TemplateGUVIntegrationStepper<Backend>* TemplateGUIntegrationDriver<Backend>
  ::GetStepper() const
{
  return fpStepper;
}

template <class Backend>
inline
TemplateGUVIntegrationStepper<Backend>* TemplateGUIntegrationDriver<Backend>
  ::GetStepper() 
{
  return fpStepper;
}

template <class Backend>
inline
int TemplateGUIntegrationDriver<Backend>
  ::GetMaxNoSteps() const
{
  return fMaxNoSteps;
}

template <class Backend>
inline
void TemplateGUIntegrationDriver<Backend>
  ::SetMaxNoSteps(int val)
{
  fMaxNoSteps= val;
}

template <class Backend>
inline
void TemplateGUIntegrationDriver<Backend>
  ::GetDerivatives(const TemplateGUFieldTrack<Backend> &y_curr, // const, INput
                         typename Backend::precision_v  charge, 
                         typename Backend::precision_v  dydx[])  // OUTput
{ 
  typename Backend::precision_v  tmpValArr[TemplateGUFieldTrack<Backend>::ncompSVEC];
  y_curr.DumpToArray( tmpValArr  );
  fpStepper -> RightHandSideVIS( tmpValArr , charge, dydx );
}


template <class Backend>
const int  TemplateGUIntegrationDriver<Backend>::fMaxStepBase = 250;  // Was 5000

#ifndef G4NO_FIELD_STATISTICS
#define GVFLD_STATS  1
#endif

#define GUDEBUG_FIELD 1

#ifdef GVFLD_STATS
#include "TH1.h"
TH1F* gHistStepsLin=0;
TH1F* gHistStepsLog=0;
TH1F* gHistStepsInit=0;
#endif

// To add much printing for debugging purposes, uncomment the following
// and set verbose level to 1 or higher value !
// #define  GUDEBUG_FIELD 1    

// ---------------------------------------------------------

//  Constructor
//
template <class Backend>
TemplateGUIntegrationDriver<Backend>
  ::TemplateGUIntegrationDriver( double  hminimum, 
                                 TemplateGUVIntegrationStepper<Backend> *pStepper,
                                 int     numComponents,
                                 int     statisticsVerbose)
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
     fOneGoodStepCallsVector(0), fOneStepCallsScalar(0),
     fVerboseLevel(3)
{  
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
  assert( pStepper != nullptr );
  
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

#ifdef GVFLD_STATS
  if( !gHistStepsLin )
       gHistStepsLin = new TH1F("hSteps", "Step size in Int-Driver", 100, 0, 100.0);
  if( !gHistStepsInit )
       gHistStepsInit = new TH1F("hSteps", "Input Step size in Int-Driver", 100, 0, 100.0);  
  if( !gHistStepsLog )
       gHistStepsLog = new TH1F("hSteps", "Log of Step size in Int-Driver", 40, -10., +10.0);
#endif

  //For track insertion
  fIndex = new int[kVectorSize];

}

//  Copy Constructor - used by Clone
template <class Backend>
TemplateGUIntegrationDriver<Backend>
::TemplateGUIntegrationDriver( const TemplateGUIntegrationDriver<Backend>& right ) 
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
   const TemplateGUVIntegrationStepper<Backend> *protStepper = right.GetStepper();
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
template <class Backend>
TemplateGUIntegrationDriver<Backend>
  ::~TemplateGUIntegrationDriver()
{ 
  if( fStatisticsVerboseLevel > 1 )
  {
    PrintStatisticsReport();
  }

  delete[] fIndex;
  // delete fpScalarDriver;
  // delete fpScalarStepper;
  // delete fpStepper;
}


template <class Backend>
TemplateGUIntegrationDriver<Backend>* 
TemplateGUIntegrationDriver<Backend>
  ::Clone() const
{
   return new TemplateGUIntegrationDriver<Backend>(*this);
}

// ---------------------------------------------------------
template<>
void
TemplateGUIntegrationDriver<vecgeom::kScalar>::OneGoodStep(  double y[],        // InOut
                             const double dydx[],
                             const double charge,         
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

  int ncompSVEC = TemplateGUFieldTrack<vecgeom::kScalar>::ncompSVEC;
  double yerr[ncompSVEC], ytemp[ncompSVEC];

  // bool verbose= false;
  // if( verbose ) std::cout << "OneGoodStep called with htry= " << htry << std::endl;
  fOneStepCallsScalar++,
  
  h = htry ; // Set stepsize to the initial trial value

  double inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  double errpos_sq=0.0;    // square of displacement error
  double errmom_sq=0.0;    // square of momentum vector difference

  int iter;

  static int tot_no_trials=0;  // thread_local
  const int max_trials=100; 

  for (iter=0; iter<max_trials ;iter++)
  {
    tot_no_trials++;
    fpStepper-> StepWithErrorEstimate( y, dydx, charge, h, ytemp, yerr);
    // fStepperCalls++;

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
       errmom_sq = sumerr_sq / magmom_sq; 
    } else {
       std::cerr << "** GUIntegrationDriver: found case of zero momentum." 
                 << " iteration=  " << iter << " h= " << h << std::endl; 
       errmom_sq = sumerr_sq;
    }
    errmom_sq *= inv_eps_vel_sq;
    errmax_sq = std::max( errpos_sq, errmom_sq ); // Square of maximum error
    
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

  if (errmax_sq > fErrcon*fErrcon)
  { 
    hnext = GetSafety()*h*std::pow(errmax_sq, 0.5*GetPowerGrow());
  }
  else
  {
    hnext = fMaxSteppingIncrease*h ; // No more than a factor of 5 increase
  }
  x += (hdid = h);

  for(int k=0;k<fNoIntegrationVariables;k++) { y[k] = ytemp[k]; }

  // std::cout<<"OneGoodStep done at iter "<<iter << " with hdid: " << hdid << std::endl;
  // std::cout<<"Resultant y is: " << y[0] << " " << y[1] << std::endl;  


  return;
}   // end of  OneGoodStep .............................


template </*class Backend*/>
bool 
TemplateGUIntegrationDriver<vecgeom::kScalar>
  ::AccurateAdvance(const TemplateGUFieldTrack<vecgeom::kScalar>& yInput,
                    const double  charge,                    
                          double  hstep,
                          double  epsilon,
                          TemplateGUFieldTrack<vecgeom::kScalar>& yOutput )
                          // typename Backend::precision_v  hinitial)
{
  // Driver for Runge-Kutta integration with adaptive stepsize control.
  // Integrate starting 'vector' y_current, over length 'hstep'
  // maintaining integration error so that relative accuracy is better
  // than 'epsilon'.
  // NOTE: The number of trial steps is limited by 'fMaxNoSteps'. Integration will 
  //       stop if this maximum is reached, and the return value will be false.
  // On return
  //  - 'yOutput' provides the values at the end of the integration interval;
  //  - the return value is 'true' if integration succeeded to the end of the interval,
  //    and 'false' otherwise.

  // using ThreeVector_v = vecgeom::Vector3D<double>;

  // std::cout<<"\n " << std::endl;
#ifdef COLLECT_STATISTICS
  constexpr double perMillion  = 1.0e-6;
  constexpr double perThousand = 1.0e-3;
  int no_warnings = 0;
#endif 

  int nstp, i;
  double x, hnext, hdid, h;

  int ncompSVEC = TemplateGUFieldTrack<vecgeom::kScalar>::ncompSVEC; //12, //to be derived from TemplateGUFieldTrack

#ifdef GUDEBUG_FIELD
  // static int dbg=1;
  // static int nStpPr=50;   // For debug printing of long integrations
  // double ySubStepStart[ncompSVEC];

  // std::cout << " AccurateAdvance called with hstep= " << hstep << std::endl;
#endif

  double  y     [ncompSVEC], 
          dydx  [ncompSVEC];
  double  ystart[ncompSVEC], 
          yEnd  [ncompSVEC]; 
  double  x1, x2;
  Bool_v succeeded(true), lastStepSucceeded;

  double startCurveLength;

  // int  noFullIntegr=0, noSmallIntegr = 0 ;
  // ThreadLocal
  // int  noGoodSteps =0 ;  // Bad = chord > curve-len 
  const  int  nvar= fNoVars;

  //  Ensure that hstep > 0
  if( hstep <= 0.0 )
  { 
    if(hstep==0.0)
    {
      std::cerr << "Proposed step is zero; hstep = " << hstep << " !" << std::endl;
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
  x2= x1 + hstep; // = x + hstep

  h = hstep; // = x2 -x1 ; or x2 - x

  x = x1;


  for (i=0; i<fNoVars; i++)  { y[i] = ystart[i]; }

  bool lastStep = false;
  nstp=1;

  // double StartPosAr[3];
  // double charge(-1.);

  while ( ((nstp++)<=fMaxNoSteps) && (x < x2) && (!lastStep) )
  {
    // StartPosAr[0]= y[0]; StartPosAr[1]= y[1];  StartPosAr[2]= y[2]; 
    
    fpStepper->RightHandSideVIS( y, charge, dydx );   // TODO: change to inline

#ifdef COLLECT_STATISTICS
    fNoTotalSteps++;
#endif 

    // Perform the Integration
    if( h > fMinimumStep )
    { 
      OneGoodStep( y, dydx, charge, x, h, epsilon, hdid, hnext) ;

      //--------------------------------------
      lastStepSucceeded= (hdid == h);   

    }

#ifdef COLLECT_STATISTICS
    if (lastStepSucceeded)  { noFullIntegr++; }
    else                    { noSmallIntegr++; }
#endif 


    // Check the endpoint
    //   ThreeVector_v EndPos( y[0], y[1], y[2] );
    // const double edx= y[0] - StartPosAr[0];
    // const double edy= y[1] - StartPosAr[1];
    // const double edz= y[2] - StartPosAr[2];
    // double endPointDist2= vecgeom::Sqrt(edx*edx+edy*edy+edz*edz) ; // (EndPos-StartPos).Mag(); 

    //Ananya: discuss. What exactly is happening here?
    bool avoidNumerousSmallSteps = (h < epsilon * hstep) || (h < fSmallestFraction * startCurveLength);
    if (avoidNumerousSmallSteps ==1)
    {
      lastStep = true;
    }
    // vecgeom::MaskedAssign(avoidNumerousSmallSteps, true, &lastStep);

    //For rest, check the proposed next stepsize 
    h = vecgeom::Max(hnext, fMinimumStep);

    // Ensure that the next step does not overshoot
    if ((x+h) > x2 )
    {
      h = x2 -x ;
    }

    if (h==0)
    {
      lastStep = true;
    }
  } 
     // Have we reached the end ?
     // --> a better test might be x-x2 > an_epsilon

  succeeded=  (x>=x2);  // If it was a "forced" last step

  for (i=0;i<nvar;i++)  { yEnd[i] = y[i]; }

  // Put back the values.
  yOutput.LoadFromArray( yEnd, fNoIntegrationVariables );
  yOutput.SetCurveLength( x );

  return succeeded;
}  // end of AccurateAdvance ...........................

// ---------------------------------------------------------

/*template <class Backend>
void
TemplateGUIntegrationDriver<Backend>
  ::WarnSmallStepSize( double hnext,
                       double hstep, 
                       double h,
                       double xDone,
                       int    nstp)
{
  static int noWarningsIssued =0;   // thread_local
  const  int maxNoWarnings =  10;   // Number of verbose warnings
  // std::ostringstream message;
  // typedef std::cerr message;
  std::cerr << " WARNING from TemplateGUIntegrationDriver::WarnSmallStepSize():  " ; // << std::endl;
  if( (noWarningsIssued < maxNoWarnings) || fVerboseLevel > 10 )
  {
    std::cerr << "The stepsize for the next iteration, " << hnext
            << ", is too small - in Step number " << nstp << "." << std::endl
              << "The minimum for the driver is " << GetHmin() << " "  // << std::endl
              << "Requested integr. length was " << hstep << " ." << std::endl
              << "The size of this sub-step was " << h     << " ." // << std::endl
              << " Integration has already done length= " << xDone << std::endl;
  }
  else
  {
    std::cerr << "Too small 'next' step " << hnext
            << ", step-no: " << nstp << std::endl
            << ", this sub-step: " << h     
            << ",  req_tot_len: " << hstep 
            << ", done: " << xDone << ", min: " << GetHmin();
  }
  // G4Exception("TemplateGUIntegrationDriver<Backend>::WarnSmallStepSize()", "GeomField1001",
  //             JustWarning, message);
  noWarningsIssued++;
}

// ---------------------------------------------------------

template <class Backend>
void
TemplateGUIntegrationDriver<Backend>
  ::WarnTooManySteps( double x1start, 
                      double x2end, 
                      double xCurrent)
{

   std::cerr << "WARNING from TemplateGUIntegrationDriver::WarnTooManySteps()" << std::endl;
   std::cerr << "The number of steps used in the Integration driver"
             << " (Runge-Kutta) is too many." << std::endl
             << "Integration of the interval was not completed !" << std::endl;

   unsigned int oldPrec= std::cerr.precision(16);
   
   std::cerr << "Only a " << (xCurrent-x1start)*100.0/(x2end-x1start)
             << " % fraction of it was done.";
   // std::cerr.setf (std::ios_base::scientific);
   // std::cerr.precision(4);
   std::cerr << "Remaining fraction= " << (x2end-xCurrent)*100.0/(x2end-x1start)      
             << std::endl;
   // std::cerr.unsetf (std::ios_base::scientific);
   std::cerr.precision(oldPrec);
   // G4Exception("TemplateGUIntegrationDriver::WarnTooManySteps()", "GeomField1001",
   //             JustWarning, message);
}

// ---------------------------------------------------------

template <class Backend>
void
TemplateGUIntegrationDriver<Backend>
  ::WarnEndPointTooFar( double endPointDist, 
                        double h , 
                        double epsilon,
                        int    dbg         )
{
  static  double maxRelError=0.0; // thread_local
  bool isNewMax, prNewMax;

  isNewMax = endPointDist > (1.0 + maxRelError) * h;
  prNewMax = endPointDist > (1.0 + 1.05 * maxRelError) * h;
  if( isNewMax ) { maxRelError= endPointDist / h - 1.0; }

  if( dbg && (h > fSurfaceTolerance) 
          && ( (dbg>1) || prNewMax || (endPointDist >= h*(1.+epsilon) ) ) )
  { 
    static int noWarnings = 0;  // thread_local
    // std::ostringstream message;
    std::cerr << "WARNING in TemplateGUIntegrationDriver::WarnEndPointTooFar()" << std::endl;
    if( (noWarnings ++ < 10) || (dbg>2) )
    {
      std::cerr << "The integration produced an end-point which "
              << "is further from the start-point than the curve length."
              << std::endl;
    }
    std::cerr << "  Distance of endpoints = " << endPointDist
              << ", curve length = " << h << std::endl
              << "  Difference (curveLen-endpDist)= " << (h - endPointDist)
              << ", relative = " << (h-endPointDist) / h 
              << ", epsilon =  " << epsilon << std::endl;
    // G4Exception("TemplateGUIntegrationDriver::WarnEndPointTooFar()", "GeomField1001",
    //             JustWarning, message);
  }
}*/

// ---------------------------------------------------------

template </*class Backend*/>
void
TemplateGUIntegrationDriver<vecgeom::kVc>
  ::OneGoodStep(       typename vecgeom::kVc::precision_v  y[],        // InOut
                 const typename vecgeom::kVc::precision_v  dydx[],
                 const typename vecgeom::kVc::precision_v  charge,
                       typename vecgeom::kVc::precision_v& x,         // InOut
                       typename vecgeom::kVc::precision_v  htry,
                       typename vecgeom::kVc::precision_v  eps_rel_max,
                       typename vecgeom::kVc::precision_v& hdid,      // Out
                       typename vecgeom::kVc::precision_v& hnext      )    // Out

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
  Double_v errmax_sq;
  Double_v h, htemp, xnew ;

  Double_v yerr [TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC], 
           ytemp[TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC];

  // std::cout << "OneGoodStep called with htry= " << htry << std::endl;
  fOneGoodStepCallsVector++;
  
  h = htry ; // Set stepsize to the initial trial value

  Double_v inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  Double_v errpos_sq=0.0;    // square of displacement error
  Double_v errmom_sq=0.0;    // square of momentum vector difference

  int iter;

  static int tot_no_trials=0;  // thread_local
  const int max_trials=100; 

  for (iter=0; iter<max_trials ;iter++)
  {
    if (true)
    {
      tot_no_trials++;
      fpStepper-> StepWithErrorEstimate(y, dydx, charge, h, ytemp, yerr);

      Double_v eps_pos = eps_rel_max * vecgeom::Max(h, fMinimumStep);  // Uses remaining step 'h'
      Double_v inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos);

      // Evaluate accuracy
      errpos_sq =  yerr[0]*yerr[0] + yerr[1]*yerr[1] + yerr[2]*yerr[2] ; 
      errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

      // Accuracy for momentum
      Double_v magmom_sq=  y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
      Double_v sumerr_sq=  yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5]; 

      vecgeom::CondAssign(magmom_sq > 0.0, sumerr_sq/magmom_sq, sumerr_sq, &errmom_sq);

      errmom_sq *= inv_eps_vel_sq;
      errmax_sq = vecgeom::Max( errpos_sq, errmom_sq ); // Square of maximum error

      // Ananya : how to break in MaskedAssign? 
      // Also, now will need to break for only one track out of 4.
      // maybe make a done? and then use that Done for every other MaskedAssign etc.. .. 
      // or wherever final value assignment is being done

      // Ananya : Need to change things here too.... As of now commenting the line below, but 
      // it is not appropriate
      // if ( errmax_sq <= 1.0 )  { break; } // Step succeeded. 

      // Step failed; compute the size of retrial Step.
      // Ananya : adding a statement. Later check the sanity or work around
      Double_v errPower = Vc::exp( (0.5*fPowerShrink)*vecgeom::Log(errmax_sq) );
      htemp = fSafetyFactor *h* errPower;

      // htemp = fSafetyFactor *h* vecgeom::Pow( errmax_sq, 0.5*fPowerShrink );

      //Can also use a max here
      h = vecgeom::Max(htemp, 0.1*h);
      // vecgeom::CondAssign(htemp >= 0.1*h, htemp, 0.1*h, &h);
  /*    if (htemp >= 0.1*h)  { h = htemp; }  // Truncation error too large,
      else  { h = 0.1*h; }                 // reduce stepsize, but no more
                                           // than a factor of 10*/
      
      xnew = x + h;

      //Ananya: Seems like using even Done won't resolve it. Can use IsFull to break 
      //if all particles happen to be stuck or something 
      //Confirm though
      if(vecgeom::IsFull(xnew==x)) {
        std::cerr << "GVIntegratorDriver::OneGoodStep:" << std::endl
         << "  Stepsize underflow in Stepper " << std::endl ;
        std::cerr << "  Step's start x=" << x << " and end x= " << xnew 
               << " are equal !! " << std::endl
               <<"  Due to step-size= " << h 
               << " . Note that input step was " << htry << std::endl;
        break;
      }
  /*    if(xnew == x)
      {
        break;
      }*/     
    }

  }
  std::cout << "GUIntDrv: 1-good-step - Loop done at iter = " << iter << std::endl;

#ifdef GVFLD_STATS
  // Sum of squares of position error // and momentum dir (underestimated)
  fSumH_lg += h; 
  fDyerrPos_lgTot += errpos_sq;
  fDyerrVel_lgTot += errmom_sq * h * h; 
#endif

  // Compute size of next Step
  Double_v errPower = Vc::exp( (0.5*GetPowerGrow())* vecgeom::Log(errmax_sq));
  hnext = GetSafety()*errPower;
  // hnext = GetSafety()*vecgeom::Pow(errmax_sq, 0.5*GetPowerGrow());
  vecgeom::MaskedAssign(errmax_sq <= fErrcon*fErrcon, fMaxSteppingIncrease*h, &hnext); // No more than a factor of 5 increase

  x += (hdid = h);

  for(int k=0;k<fNoIntegrationVariables;k++) { y[k] = ytemp[k]; }

  return;
}   // end of  OneGoodStep .............................

//----------------------------------------------------------------------

#define SQR(a)   ((a)*(a))

// QuickAdvance just tries one Step - it does not ensure accuracy
template <class Backend>//
typename Backend::bool_v  
TemplateGUIntegrationDriver<Backend>
  ::QuickAdvance( TemplateGUFieldTrack<Backend>&       y_posvel,    // InOut
                  const typename Backend::precision_v  dydx[],
                  const typename Backend::precision_v  charge,      // In
                        typename Backend::precision_v  hstep,       // In
#ifdef USE_DCHORD
                        typename Backend::precision_v& dchord_step,
#endif
                        typename Backend::precision_v& dyerr_pos_sq,
                        typename Backend::precision_v& dyerr_mom_rel_sq )  
{
  typedef typename Backend::precision_v Double_v;
//  typedef typename Backend::bool_v      Bool_v;
  // Double_v dyerr_pos_sq, dyerr_mom_rel_sq;  
  Double_v yerr_vec[TemplateGUFieldTrack<Backend>::ncompSVEC],
           yarrin  [TemplateGUFieldTrack<Backend>::ncompSVEC], 
           yarrout [TemplateGUFieldTrack<Backend>::ncompSVEC]; 
  Double_v s_start;
  Double_v dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

  static int no_call=0;  // thread_local 
  no_call ++; 

  // Move data into array
  y_posvel.DumpToArray( yarrin );      //  yarrin  <== y_posvel 
  s_start = y_posvel.GetCurveLength();

  // Do an Integration Step
  fpStepper-> StepWithErrorEstimate( yarrin, dydx, charge, hstep, yarrout, yerr_vec ) ; 
  //          *********************

#ifdef USE_DCHORD  
  // Estimate curve-chord distance
  dchord_step= fpStepper-> DistChord();
  //                       *********
#endif

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

#ifdef RETURN_A_NEW_STEP_LENGTH
  // The following step cannot be done here because "eps" is not known.
  dyerr_len = vecgeom::Sqrt( dyerr_len_sq ); 
  dyerr_len_sq /= epsilon ;

  // Set suggested new step
  hstep= ComputeNewStepSize( dyerr_len, hstep);
#endif

  return true;
}

// --------------------------------------------------------------------------
#ifdef QUICK_ADV_ARRAY_IN_AND_OUT
template <class Backend>
typename Backend::bool_v  
TemplateGUIntegrationDriver<Backend>
  ::QuickAdvance(       typename Backend::precision_v     yarrin[],    // In
                  const typename Backend::precision_v     dydx[],
                  const typename Backend::precision_v     charge,      // In
                        typename Backend::precision_v     hstep,       // In
                        typename Backend::precision_v     yarrout[],
                        typename Backend::precision_v&    dchord_step,
                        typename Backend::precision_v&    dyerr )      // In length
{
   std::cerr << "ERROR in TemplateGUIntegrationDriver::QuickAdvance()" << std::endl;
   std::cerr << "      Method is not yet implemented." << std::endl;

   //            FatalException, "Not yet implemented.");
   dyerr = dchord_step = hstep * yarrin[0] * dydx[0];
   yarrout[0]= yarrin[0];
   exit(1); 
}
#endif 

// --------------------------------------------------------------------------

//  This method computes new step sizes - but does not limit changes to
//  within  certain factors
// 
template <class Backend>
typename Backend::precision_v 
TemplateGUIntegrationDriver<Backend>
  ::ComputeNewStepSize( typename Backend::precision_v  errMaxNorm,    // max error  (normalised)
                        typename Backend::precision_v  hstepCurrent)  // current step size
{
  typename Backend::precision_v hnew;

  // Compute size of next Step for a failed step
  hnew = GetSafety()*hstepCurrent*vecgeom::Pow(errMaxNorm,GetPowerShrink()) ;
  
  vecgeom::MaskedAssign(errMaxNorm <= 1.0 /*&& errMaxNorm >0.0*/, GetSafety()*hstepCurrent*vecgeom::Pow(errMaxNorm,GetPowerGrow()), &hnew );

  vecgeom::MaskedAssign(errMaxNorm <= 0.0, fMaxSteppingIncrease * hstepCurrent, &hnew);

  hnew = fMaxSteppingIncrease * hstepCurrent; 


/*  if(errMaxNorm > 1.0 )
  {
    hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPowerShrink()) ;
  } 
  else if(errMaxNorm > 0.0 ) {
    hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPowerGrow()) ;
  } 
  else {
    hnew = fMaxSteppingIncrease * hstepCurrent; 
  }*/

  return hnew;
}

// ---------------------------------------------------------------------------

// This method computes new step sizes limiting changes within certain factors
// 
// It shares its logic with AccurateAdvance.
// They are kept separate currently for optimisation.
//
template <class Backend>
typename Backend::precision_v 
TemplateGUIntegrationDriver<Backend>
  ::ComputeNewStepSize_WithinLimits( typename Backend::precision_v  errMaxNorm,    // max error  (normalised)
                                     typename Backend::precision_v  hstepCurrent)  // current step size
{
  typename Backend::precision_v hnew, htemp;

  hnew  = fMaxSteppingIncrease * hstepCurrent;

  htemp = GetSafety()*hstepCurrent*vecgeom::Pow(errMaxNorm,GetPowerShrink()) ; 
  //Size of next step for failed step
  vecgeom::MaskedAssign(errMaxNorm>1.0 && htemp > hnew, htemp, &hnew);

  htemp= GetSafety()*hstepCurrent*vecgeom::Pow(errMaxNorm,GetPowerGrow());
  //Size of next step for successful step
  vecgeom::MaskedAssign(errMaxNorm<=1.0 && errMaxNorm > fErrcon, htemp, &hnew);


/*
  // Compute size of next Step for a failed step
  if (errMaxNorm > 1.0 )
  {
    // Step failed; compute the size of retrial Step.
    hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPowerShrink()) ;
  
    if (hnew < fMaxSteppingDecrease*hstepCurrent)
    {
      hnew = fMaxSteppingDecrease*hstepCurrent ;
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
     { hnew = fMaxSteppingIncrease * hstepCurrent; }
  }*/

  return hnew;
}

// ---------------------------------------------------------------------------
/*template <class Backend>
void 
TemplateGUIntegrationDriver<Backend>
::PrintStatus( const double*   StartArr,  
                     double    xstart,
               const double*   CurrentArr, 
                     double    xcurrent,
                     double    requestStep, 
                     int       subStepNo  )
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
template <class Backend>
void TemplateGUIntegrationDriver<Backend>
  ::PrintStatus( const TemplateGUFieldTrack<Backend>&  StartFT,
                 const TemplateGUFieldTrack<Backend>&  CurrentFT, 
                 double             requestStep, 
                 int                subStepNo  )
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
              << " TemplateGUIntegrationDriver: Current Position  and  Direction" << " "
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
      
    }

    if( verboseLevel <= 3 )
    {
      std::cout.precision(noPrecision);
      PrintStat_Aux( CurrentFT, requestStep, step_len, 
                     subStepNo, subStepSize, DotStartCurrentVeloc );
      
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
template <class Backend>
void TemplateGUIntegrationDriver<Backend>
  ::PrintStat_Aux( const TemplateGUFieldTrack<Backend>&  aGUFieldTrack,
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
***/

// ---------------------------------------------------------------------------
template <class Backend>
void 
TemplateGUIntegrationDriver<Backend>
  ::PrintStatisticsReport()
{
  int noPrecBig= 6;
  int oldPrec= std::cout.precision(noPrecBig);

  std::cout << "TemplateGUIntegrationDriver Statistics of steps undertaken. " << std::endl;
  std::cout << "TemplateGUIntegrationDriver<Backend>: Number of Steps: "
         << " Total= " <<  fNoTotalSteps
         << " Bad= "   <<  fNoBadSteps 
         << " Small= " <<  fNoSmallSteps 
         << " Non-initial small= " << (fNoSmallSteps-fNoInitialSmallSteps)
         << std::endl;

  std::cout << " Calles to Step methods:  1-Good/Vector " << fOneGoodStepCallsVector <<  " 1-step/Scalar: " << fOneStepCallsScalar << std::endl;
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
template <class Backend>
void TemplateGUIntegrationDriver<Backend>
  ::SetSmallestFraction(double newFraction)
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



// #ifdef NEWACCURATEADVANCE
template <class Backend>
void 
TemplateGUIntegrationDriver<Backend>
  ::SetNTracks(int nTracks)
// Set fNTracks 
{
  fNTracks = nTracks;
}

// #define DEBUG

template </*class Backend*/>
void 
TemplateGUIntegrationDriver<vecgeom::kVc>
  ::InitializeLanes( /*const*/              FieldTrack    yInput[],
                       const                double        hstep [],
                       const                double        chargeIn[],
                     typename vecgeom::kVc::precision_v   yLane[],
                     typename vecgeom::kVc::precision_v & hStepLane,
                     typename vecgeom::kVc::precision_v & chargeLane,                     
                     typename vecgeom::kVc::precision_v & startCurveLength)
// Initialization step for AccurateAdvance/
// Converts input scalar stream to acceptable form of Vc vectors
// for vector processing in OneStep
{
#ifdef DEBUG
  std::cout<<"----Initializing lanes for AccurateAdvance----"<<std::endl;
#endif

  double yStartScalar[fNoVars];   // fNoVars or fNoIntegrationVariables? Ask : Ananya
  for (int j = 0; j < kVectorSize; ++j)
  {
    fIndex          [j] = j;
    hStepLane       [j] = hstep[j];
    chargeLane      [j] = chargeIn[j];    
    startCurveLength[j] = yInput[j].GetCurveLength();

    yInput[j].DumpToArray(yStartScalar);

    for (int i = 0; i < fNoVars; ++i)
    {
      yLane[i][j] = yStartScalar[i]; 
    }
  }
} // End of InitializeLanes function

// InsertNewTrack:  Fill working arrays ('lanes') with next input tracks, for finished (succeeded) lanes
//  @param   yInput    :  FieldTrack[]    - containter of input tracks
//  @param   hstep     :  const double[]  - total integration length requested                      
//  @param   chargeIn  :  const double[]  - charge of input tracks
//  @param   currIndex :  int             - lane to refill
//  @param   trackNextInput: int&         - counter to next available location of input (updated)
//  @param   succeeded :  bool[]          - flag of finished lanes
//  @param   y         :  Double_v[N]     - Working 'arrays' of integrands  M_vector_len x  N_components
//  @param   hStepLane :  Double_v&       - Working value, total length needed (per lane)
//  @param   chargeLane:  Double_v&       - Working value, charge (per lane)
//  @param   startCurveLength:  Double_v& - Working value, initial curve length (per lane)

template </*class Backend*/>
bool 
TemplateGUIntegrationDriver<vecgeom::kVc>
  ::InsertNewTrack( /*const*/ FieldTrack                    yInput[],
                        const double                        hstep[],
                        const double                        chargeIn[],                    
                              int                           currIndex,
                              int&                          trackNextInput,
                              bool                          succeeded[],
                        typename vecgeom::kVc::precision_v  y[],
                        typename vecgeom::kVc::precision_v &hStepLane,
                        typename vecgeom::kVc::precision_v &chargeLane,                    
                        typename vecgeom::kVc::precision_v &startCurveLength )
// Inserts a new track whenever a lane is finished. 
// returns isDoneLane = true for h<=0 case, false otherwise
// because in former case, no further work is required
{
#ifdef DEBUG
  std::cout<<"----Inserting New Track "<< trackNextInput << " at position "<< currIndex <<std::endl;
#endif 

  bool found = false; // to get the while loop starting
  while(trackNextInput < fNTracks && !found)
  {

    // Ensure that hstep > 0
    double hStepNext = hstep[trackNextInput];
    if (hStepNext <= 0)
    {
      if (hStepNext == 0 )
      {
        std::cerr << "Proposed step is zero; hstep = " << hStepNext << " !" << std::endl;
        // Success step
        // succeeded is initialized with true, hence no assignment needed here
      }
      else
      {
        std::cerr << "Invalid run condition." << std::endl
                  << "Proposed step is negative; hstep = " << hStepNext << "." << std::endl;
        succeeded[trackNextInput] = false; //the final bool array to be returned 
      }
    }

    else
    {
      found = true;
      double yScalar[fNoVars];
      yInput[trackNextInput].DumpToArray(yScalar);
      for (int i = 0; i < fNoVars; ++i)
      {
        y[i][currIndex] = yScalar[i];    // y[0] = Double_v of '0'/x coordinates
      }
      fIndex          [currIndex] = trackNextInput;
      hStepLane       [currIndex] = hstep[trackNextInput];
      chargeLane      [currIndex] = chargeIn[trackNextInput];      
      startCurveLength[currIndex] = yInput[trackNextInput].GetCurveLength();
    }

    trackNextInput++;
  }
  
  return found;

} // End of InsertNewTrack function

template </*class Backend*/>
void 
TemplateGUIntegrationDriver<vecgeom::kVc>
  ::StoreOutput(const typename vecgeom::kVc::precision_v y[],
                const typename vecgeom::kVc::precision_v x,
                      FieldTrack                    yOutput[],
                      int                           currIndex,
                      double                        hstep[],
                      bool                          succeeded[] )
// Called whenever a lane is finished.
// Stores value of succeeded in the bool[nTracks]
// Stores final curve length and end position and momentum
// in yOutput (scalar form)
// currIndex is the index of finished lane in Vc vector 
// hstep argument given because no storage needed for 
// h<=0 case 
{
#ifdef DEBUG
  cout<<"----Storing Output at position: "<< currIndex << std::endl;
#endif 

  int    absoluteIndex = fIndex[currIndex]; // might be sent directly to StoreOutput as well
  double hStepOriginal = hstep [absoluteIndex];

  if(hStepOriginal == 0.0)
  {
    succeeded[absoluteIndex] = true;
  }
  else if (hStepOriginal < 0.0)
  {
    succeeded[absoluteIndex] = false;
  }
  else
  {
    // need to get a yEnd : scalar array
    double yEnd[fNoVars]; // Confirm size //ncompSVEC? : Ananya 
    for (int i = 0; i < fNoIntegrationVariables; ++i)
    {
      yEnd[i] = y[i][currIndex]; // Constant col no., varying row no. for required traversal
    }
    yOutput[absoluteIndex].LoadFromArray ( yEnd, fNoIntegrationVariables);
    yOutput[absoluteIndex].SetCurveLength( x[currIndex]); // x is a double_v variable
  }

} // End of StoreOutput function


// OneStep - integrate each lanes by one successful step (?)
//   @param   y[]     - Integrand (x, p)                - In/Out
//   @param   dydx[]  - Derivative                      - In
//   @param   charge  - particle charge                 - In (const)
//   @param   x       - independent variable  ( in= initial value, out = current value)
//   @param   htry    - proposed size for first step    - In (const)
//   @param   epsRel  - required relative accuracy      - In (const)
//   @param   hdid    - step size achieved in call      - Out    ????
//   @param   hnext   - proposed next step size         - Out

template </*class Backend*/>
void
TemplateGUIntegrationDriver<vecgeom::kVc>  // TemplateGUIntegrationDriver<Backend>
  ::OneStep(       typename vecgeom::kVc::precision_v  y[],        // InOut
             const typename vecgeom::kVc::precision_v  dydx[],
             const typename vecgeom::kVc::precision_v  charge,      // In                   
             const typename vecgeom::kVc::precision_v  eps_rel_max,
                   typename vecgeom::kVc::precision_v& x,           // InOut
                   typename vecgeom::kVc::precision_v  htry,
                   typename vecgeom::kVc::precision_v& hdid,        // Out
                   typename vecgeom::kVc::precision_v& hnext      ) // Out

// Derived from OneGoodStep
//   This version attempts to make one good step in each of the vector 'lanes'
//   If the trial step is not successful, it will reduce the step-size.
//   It will finish in a lane, if:
//     - it has succeeded for the current step-size  (original or revised)
//     - there is underflow in the trial step size 'h',
//     - the maximum number of iterations is reached ( max_trials - local const)

// Driver for one Runge-Kutta Step with monitoring of local truncation error
// to ensure accuracy and adjust stepsize. Input are dependent variable
// array y[0,...,5] and its derivative dydx[0,...,5] at the
// starting value of the independent variable x . Also input are stepsize
// to be attempted htry, and the required accuracy eps. On output y and x
// are replaced by their new values, hdid is the stepsize that was actually
// accomplished, and hnext is the estimated next stepsize. 
// This is similar to the function rkqs from the book:
// Numerical Recipes in C:p. 719

{
#ifdef PARTDEBUG
  if (partDebug) { std::cout<<"\n"<<std::endl; }
#endif 
  Double_v h, htemp, xnew, errmax_sq;
  Double_v yerr [ncompSVEC],
           ytemp[ncompSVEC];

  h = htry ; // Set stepsize to the initial trial value

  Double_v inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);
  Double_v errpos_sq=0.0;    // square of displacement error
  Double_v errmom_sq=0.0;    // square of momentum vector difference

  int iter;

  // static int tot_no_trials=0;  // thread_local
  const  int max_trials=100; 

  Bool_v   finished(false); // [kVectorSize] = {false, false, false, false};  // Care: AVX-only

  Double_v hFinal, hnextFinal, xFinal, hdidFinal, errmax_sqFinal;
  Double_v yFinal[ncompSVEC];
  Bool_v errMaxLessThanOne(false), hIsZeroCond(false);

  for (iter = 0; iter < max_trials; iter++)
  {
    if (  !vecgeom::IsFull(hIsZeroCond || errMaxLessThanOne) )
    {
      // tot_no_trials++;
      fpStepper-> StepWithErrorEstimate(y, dydx, charge, h, ytemp, yerr);
      fStepperCalls++;

#ifdef PARTDEBUG
      if (partDebug) {
        std::cout<< "----h is: " << h[0] << " at iter: " << iter << std::endl;
        std::cout<< "----yerr is [0]: " << yerr[0][0] <<" [1]: "<<yerr[1][0]
                 << " [2]: " << yerr[2][0] << " [3]: " << yerr[3][0]
                 << " [4]: " << yerr[4][0] << " [5]: " << yerr[5][0] << std::endl;
      }
#endif 
      Double_v eps_pos = eps_rel_max * vecgeom::Max(h, fMinimumStep);  // Uses remaining step 'h'
      Double_v inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos);

      // Evaluate accuracy
      errpos_sq =  yerr[0]*yerr[0] + yerr[1]*yerr[1] + yerr[2]*yerr[2] ; 
      errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

      // Accuracy for momentum
      Double_v magmom_sq=  y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
      Double_v sumerr_sq=  yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5]; 

      // errmom_sq =  (magmom_sq > 0.0 ) ? sumerr_sq/magmom_sq : sumerr_sq;
      vecgeom::CondAssign(magmom_sq > 0.0, sumerr_sq/magmom_sq, sumerr_sq, &errmom_sq);
      errmom_sq *= inv_eps_vel_sq;
      errmax_sq = vecgeom::Max( errpos_sq, errmom_sq ); // Square of maximum error
#ifdef PARTDEBUG
      if (partDebug) {
        std::cout<< "----eps_pos is: "<< eps_pos[0] << std::endl;
        std::cout<< "----inv_eps_pos_sq is: "<< inv_eps_pos_sq[0] << std::endl;
        std::cout<< "----errmom_sq is: "<< errmom_sq[0] << std::endl;
        std::cout<< "----errpos_sq is: "<< errpos_sq[0] << std::endl;
        std::cout<< "----errmax_sq is: "<< errmax_sq[0] << std::endl;
      }
#endif 
      errMaxLessThanOne = ( errmax_sq <=1.0 );

      // int numOngoing = ?
      if ( !vecgeom::IsEmpty(errMaxLessThanOne) )
      {
        //  Some lanes succeeded -- store their output
        for (int i = 0; i < kVectorSize; ++i)
        {
          // bool justFinished[i] = ( errMaxLessThanOne[i] && !finished[i] )           
          if ( errMaxLessThanOne[i] && !finished[i] )
          {
            /* StoreFinalValues() */
            // - Could use several MaskedAssigns instead             
            finished      [i] = true;
            hFinal        [i] = h[i];
            errmax_sqFinal[i] = errmax_sq[i];
            for (int j = 0; j < ncompSVEC; ++j)
            {
              yFinal[j][i] = ytemp[j][i];
            }
          }
        }
      }
      if ( vecgeom::IsFull(errMaxLessThanOne) )  { break; } // Step succeeded. 

      Double_v  hLast= h;
      // Step failed; compute the size of retrial Step.
      // Ananya : adding a statement. Later check the sanity or work around
#if 0
      Double_v errPower = Vc::exp( (0.5*fPowerShrink)*vecgeom::Log(errmax_sq) ); 
      htemp = fSafetyFactor * h * errPower;
      // htemp = fSafetyFactor *h* vecgeom::Pow( errmax_sq, 0.5*fPowerShrink );
#else
      // Use loop to avoid cost of unused power operations (which are expensive)
      for (int i = 0; i < kVectorSize; ++i)
      {
        if ( ! finished[i] )
        {
          htemp[i] = fSafetyFactor * h[i] * std::pow( errmax_sq[i], 0.5*fPowerShrink);
        }
      }
#endif
      h = vecgeom::Max(htemp, 0.1*h);
      xnew = x + h;

      // enum IntegrationStatus { kStillWorking= 0,  kDone= 1, kStepUnderflow=2, kTooManySteps= 3} finished;
      
      hIsZeroCond = (xnew == x) &&  !finished;
      if ( !vecgeom::IsEmpty(hIsZeroCond) )
      {
        for (int i = 0; i < kVectorSize; ++i)
        {
          // Probably could use several MaskedAssigns as well
          if ( hIsZeroCond[i] ==1 ) // && !finished[i] )
          {
            /* StoreFinalValues() */
            finished      [i] = true; // kUnderflow;
            hFinal        [i] = h[i];
            errmax_sqFinal[i] = errmax_sq[i];
            for (int j = 0; j < ncompSVEC; ++j)
            {
              yFinal[j][i] = ytemp[j][i];
            }

            if( 1 ) { //  ( verbose ) {
               const std::string  MethodName("TemplateGUIntegrationDriver::OneStep");
               double DyDxI[ncompSVEC];
               for (int j = 0; j < ncompSVEC; ++j) { DyDxI[j] = dydx[i][j]; }
               ReportStepUnderflow( MethodName, x[i], xnew[i], htry[i], hLast[i], h[i], errpos_sq[i], errmom_sq[i],
                                    DyDxI, ncompSVEC );
            }

          }
        }
      }              
      if( vecgeom::IsFull(finished) ) {
         break;
      }
    } // end of if (  !vecgeom::IsFull(hIsZeroCond ...
  } // end of 'for iter'


#ifdef PARTDEBUG
  if (partDebug)
  {
    std::cout << "TemplateGUIntDrv: 1--step - Loop done at iter = " << iter << " with htry= " << htry <<std::endl;
  }
#endif 

  h         = hFinal;
  errmax_sq = errmax_sqFinal;

  // Compute size of next Step
  Double_v errPower = Vc::exp( (0.5*GetPowerGrow())* vecgeom::Log(errmax_sq));
  hnext = GetSafety() * h * errPower;
  // hnext = GetSafety()*vecgeom::Pow(errmax_sq, 0.5*GetPowerGrow());
  vecgeom::MaskedAssign(errmax_sq <= fErrcon*fErrcon, fMaxSteppingIncrease*h, &hnext); // No more than a factor of 5 increase

  hdid = h;
  x += h;

  for(int k=0; k < fNoIntegrationVariables; k++) {
     y[k] = yFinal[k];
  }

#ifdef PARTDEBUG
  if (partDebug) {
    std::cout<< "OneStep done - hdid= " << h << " and hnext= " << hnext << std::endl;
  }
#endif
  return;
}   // end of  OneStep .............................


// KeepStepping - integrate all lanes until one is completely finished
//   @param   y[]     - Integrand (x, p)                - In/Out
//   @param   dydx[]  - Derivative                      - In/Out (tbc)
//   @param   x       - independent variable  ( in= initial value, out = current value)
//   @param   charge  - particle charge                 - In (const)
//   @param   htry    - proposed size for first step    - In (const)
//   @param   hTotalLane - total size of integration    - In (const)
//   @param   epsRel  - required relative accuracy      - In (const)
//   @param   hnext   - proposed next step size         - Out
//   @param   hDone   -  total step size achieved       - Out 
//   @param   stepsLane - number of steps taken in lane - In/Out
// Supressed:
//   @param   hdid    - step size achieved in call, in last step ???  - Out  
// Extension: (tbc)
//   @param   lastWork - This is the last set of lanes. Return when all is done. - In (const)

template </*class Backend*/>
void
TemplateGUIntegrationDriver<vecgeom::kVc>
  ::KeepStepping(       typename vecgeom::kVc::precision_v  y[],        // InOut
                        typename vecgeom::kVc::precision_v  dydx[],     // InOut ?
                        typename vecgeom::kVc::precision_v& x,          // InOut
                  const typename vecgeom::kVc::precision_v  charge,
                  const typename vecgeom::kVc::precision_v  htry,
                  const typename vecgeom::kVc::precision_v  hTotalLane,
                  const typename vecgeom::kVc::precision_v  epsRel,
                //      bool                                lastWork,  
                        typename vecgeom::kVc::precision_v& hnext,     // Out
                        typename vecgeom::kVc::precision_v& hDone,
                        typename vecgeom::kVc::precision_v& numStepsLane  // InOut
     )

// Build starting from OneGoodStep

//  All lanes step together.
//  When one lane is completely done, the method returns with
//   the results:
//    - final values at least one lane, (or abandon due to underflow)
//    - intermediate integration results in the remaining lanes

{
#ifdef PARTDEBUG
  if (partDebug)
  {
    std::cout<<"\n"<<std::endl;
  }
#endif 

  Double_v errmax_sq;
  Double_v h, htemp, xnew ;

  Double_v yerr [TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC], 
           ytemp[TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC];

  h = htry ; // Set stepsize to the initial trial value

  Double_v inv_eps_vel_sq = 1.0 / (epsRel*epsRel);

  Double_v errpos_sq = 0.0;    // square of displacement error
  Double_v errmom_sq = 0.0;    // square of momentum vector difference

  int iter;


  static int tot_no_trials = 0;  // thread_local
  const  int max_trials    = 100; 

  // int finished[kVectorSize] = {0}; // This makes all elements of array 0 
  // for (int i = 0; i < kVectorSize; ++i ) finished[i] = 0;

  Double_v hFinal(0.), hnextFinal, xFinal, errmax_sqFinal;
  // Double_v hdidFinal;
  Double_v yFinal[TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC]; // = y[]
  for (int i = 0; i < TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC; ++i)
  {
    yFinal[i] = y[i];
  }

  Bool_v errMaxLessThanOne(false), hIsZeroCond(false);

  Double_v x2 = x + (hTotalLane - hDone);   // Is precision loss for x2 an issue ?
  // Double_v x2 = x + htry;
  Bool_v  errMaxLessThanOneLocal(false), hIsZeroCondLocal(false);
  // int htryExhausted[kVectorSize] = {0};
  Bool_v   htryExhausted(false);
  // Double_v tempCharge(+1.);

  for (iter = 0; iter < max_trials; iter++)
  {
    if ( !vecgeom::IsFull(hIsZeroCond || errMaxLessThanOne) )

    // if ( vecgeom::IsEmpty(htryExhausted) && !vecgeom::IsFull(hIsZeroCond || errMaxLessThanOne) )
    {
    #ifdef PARTDEBUG  
      if (partDebug)
      {
        Bool_v hZeroOrErrCond = hIsZeroCond || errMaxLessThanOne;
        std::cout<< "hZeroOrErrCond is: "<< hZeroOrErrCond << std::endl;
      }

    #endif  

      tot_no_trials++;

      fpStepper-> RightHandSideVIS( yFinal, charge, dydx );
      fpStepper-> StepWithErrorEstimate( yFinal, dydx, charge, h, ytemp, yerr);
      fStepperCalls++;

      numStepsLane++;  // All lanes are progress

#ifdef DEBUG
      if (partDebug)
      {
        std::cout<< "\n----yerr is: " << yerr[0] <<" "<<yerr[1]<<" "<<yerr[2] << std::endl;
      }
#endif 
      Double_v eps_pos = epsRel * vecgeom::Max(h, fMinimumStep);  // Uses remaining step 'h'
      Double_v inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos);

      // Evaluate accuracy
      errpos_sq  =  yerr[0]*yerr[0] + yerr[1]*yerr[1] + yerr[2]*yerr[2] ; 
      errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

      // Accuracy for momentum
      Double_v magmom_sq =  y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
      Double_v sumerr_sq =  yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5]; 

      vecgeom::CondAssign(magmom_sq > 0.0, sumerr_sq/magmom_sq, sumerr_sq, &errmom_sq);


      errmom_sq *= inv_eps_vel_sq;
      errmax_sq  = vecgeom::Max( errpos_sq, errmom_sq ); // Square of maximum error
#ifdef DEBUG
      if (partDebug)
      {
        std::cout<< "----eps_pos is       : "<< eps_pos   << std::endl;
        std::cout<< "----inv_eps_pos_sq is: "<< eps_pos   << std::endl;
        std::cout<< "----errmom_sq is     : "<< errmom_sq << std::endl;
        std::cout<< "----errpos_sq is     : "<< errpos_sq << std::endl;
        std::cout<< "----errmax_sq is     : "<< errmax_sq << std::endl;
      }
#endif 
      xnew = x + h; // Ananya : adding here to give an appropriate value to xnew
                    // Discuss with john if appropriate 
      errMaxLessThanOneLocal = ( errmax_sq <= 1.0 );
      if ( !vecgeom::IsEmpty(errMaxLessThanOneLocal) )
      {
        for (int i = 0; i < kVectorSize; ++i)
        {
          // Probably could use several MaskedAssigns as well
          if ( errMaxLessThanOneLocal[i] ==1 && htryExhausted[i] == false )
          {

            //----- StoreFinalValues() ---- 
            // finished         [i] = -1;
            errMaxLessThanOne[i] = 1;
            xnew             [i] = x[i] + hFinal[i];
            if (xnew[i] < x2[i] )
            {
              // stay in loop 
              // But how? 
              // htryExhausted[i] = -1; // i.e. still in working condition
              // set it to -1 otherwise and store output 
              // but if at some point x > x2 , then? what about keeping output
              // updated at every xnew < x2 step? 
              // too much storage to and forth... more time taken in data transfer 
              hFinal         [i] += h[i];

            #ifdef PARTDEBUG
              if (partDebug)
              {
                std::cout<<"hFinal["<<i<<"] is: "<<hFinal[i]<<" (errMaxLessThanOne Loop)(iter= )"<<iter<<std::endl;
              }
            #endif 
              errmax_sqFinal [i] = errmax_sq[i];
              for (int j = 0; j < TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC; ++j)
              {
                yFinal[j][i] = ytemp[j][i];
              }
            }
            else
            {
              htryExhausted[i] = true;
            }

            if (xnew[i] <= x2[i])
            {
              hDone[i] += h[i]; // Diff. from hFinal because can be carried over 
                                          // from prev. KeepStepping. 
            }
          }
        }
      }
    #ifdef PARTDEBUG
      if ( vecgeom::IsFull(errMaxLessThanOneLocal) )  
      {
        if (partDebug)
        {
          std::cout<<"errMaxLessThanOneLocal is full. "<<std::endl; 
        } 
      } // Step succeeded. 

    #endif   
      // if ( vecgeom::IsFull(errMaxLessThanOneLocal) )  { break; } // Step succeeded. 

      // Step failed; compute the size of retrial Step.
      Double_v errPower = Vc::exp( (0.5*fPowerShrink)*vecgeom::Log(errmax_sq) ); 
      htemp = GetSafety() *h* errPower;


      h = vecgeom::Max(htemp, 0.1*h);

      h = vecgeom::Min(h, hTotalLane - hDone);
      
      xnew = x + h;

      hIsZeroCondLocal = (xnew == x);
      if ( !vecgeom::IsEmpty(hIsZeroCondLocal) )
      {
        for (int i = 0; i < kVectorSize; ++i)
        {
          // Probably could use several MaskedAssigns as well
          if ( hIsZeroCondLocal[i] ==1 && htryExhausted[i] == false )
          {

            /* StoreFinalValues() */
            // finished   [i] = -1;
            hIsZeroCond[i] = 1;
            xnew       [i] = x[i] + hFinal[i];
            if (xnew[i] < x2[i] )
            {
              // stay in loop 
              // But how? 
              // htryExhausted[i] = -1; // i.e. still in working condition
              // set it to -1 otherwise and store output 
              // but if at some point x > x2 , then? what about keeping output
              // updated at every xnew < x2 step? 
              // too much storage to and forth... more time taken in data transfer 
              hFinal   [i] += h[i];

            #ifdef PARTDEBUG
              if (partDebug)
              {
                std::cout<<"hFinal["<<i<<"] is: "<<hFinal[i]<<" (hIsZero Loop)(iter= )"<<iter<<std::endl;
              }
            #endif 
              errmax_sqFinal [i] = errmax_sq[i];
              for (int j = 0; j < TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC; ++j)
              {
                yFinal[j][i] = ytemp[j][i];
              }
            }
            else
            {
              htryExhausted[i] = true;
            }

            if (xnew[i] <= x2[i])
            {
              hDone[i] += h[i]; // Diff. from hFinal because can be carried over 
                                          // from prev. KeepStepping. 
            }
          }
        }
      }


/*      if ( !vecgeom::IsEmpty(hIsZeroCondLocal) )
      {
        for (int i = 0; i < kVectorSize; ++i)
        {
          if ( hIsZeroCondLocal[i] ==1 && htryExhausted[i] != true )
          {
            //---- StoreFinalValues() ----
            finished      [i]  = -1;
            hIsZeroCond   [i]  = true;

            if (xnew[i]<= x2[i])
            {
              hFinal        [i] += h[i];
            #ifdef PARTDEBUG
              std::cout<<"hFinal["<<i<<"] is: "<<hFinal[i]<<" (hIsZero Loop)"<<std::endl;
            #endif 
              errmax_sqFinal[i]  = errmax_sq[i];
              for (int j = 0; j < TemplateGUFieldTrack<Backend>::ncompSVEC; ++j)
              {
                yFinal[j][i] = ytemp[j][i];
              }             
            }
            else
            {
              htryExhausted[i] = false;
            }

          }
        }
      }  */

      if(vecgeom::IsFull(xnew==x)) {
        std::cerr << "GVIntegratorDriver::OneStep:" << std::endl
         << "  Stepsize underflow in Stepper " << std::endl ;
        std::cerr << "  Step's start x=" << x << " and end x= " << xnew 
               << " are equal !! " << std::endl
               <<"  Due to step-size= " << h 
               << " . Note that input step was " << htry << std::endl;
        break;
      }   
    }
  }

  h         = hFinal;
  errmax_sq = errmax_sqFinal;

  // Compute size of next Step
  Double_v errPower = Vc::exp( (0.5*GetPowerGrow())* vecgeom::Log(errmax_sq));
  hnext = GetSafety() * h * errPower;
  vecgeom::MaskedAssign(errmax_sq <= fErrcon*fErrcon, fMaxSteppingIncrease*h, &hnext); // No more than a factor of 5 increase

  // std::cout<<"fPowerShrink is: "<<1/fPowerShrink<<" and fPowerGrow is: "<<1/GetPowerGrow()<< std::endl;
  // hdid = h;
  x += h;

  for(int k=0;k<fNoIntegrationVariables;k++) { y[k] = yFinal[k]; }

#ifdef PARTDEBUG
  if (partDebug)
  {
    std::cout << "TemplateGUIntDrv: 1--step - Loop done at iter = " << iter << " with htry= " << htry <<std::endl;
    std::cout<< " hdid= " << h << " and " << "hnext= " << hnext <<  std::endl;
    std::cout<< "htryExhausted is: "<< htryExhausted << std::endl;
  }
#endif 

  return;
}   // end of  KeepStepping .............................

// Advance integration of multiple tracks, ensuring small-enough error
//  @param  yInput[]   : input state - position, momentum (integrands) for each track
//  @param  chargeIn[] : charge for each track
//  @param  hstep[]    : requested step size
//  @param  epsilon    : relative tolerance for integration
//  @param  yOutput[]  : output state for each track                           - OUT
//  @param  nTracks    : number of tracks
//  @param  succeeded  : was end of integration interval reach for this track  - OUT
//

#ifdef NEWACCURATEADVANCE
template </*class Backend*/>
void
TemplateGUIntegrationDriver<vecgeom::kVc>
  ::AccurateAdvance(/*const*/ FieldTrack yInput[],
                      const   double     chargeIn[],
                              double     hStepIn[],
                              double     epsilon,
                              FieldTrack yOutput[],
                              int        nTracks,
                              bool       succeeded[] )
{
  // Built on original AccurateAdvance. Takes buffer stream of nTracks
  // Converts them to Vc vectors for processing
  // Inserts new track when processing for a lane is finished.

  // Driver for Runge-Kutta integration with adaptive stepsize control.
  // Integrate starting 'vector' y_current, over length 'hStepIn'
  // maintaining integration error so that relative accuracy is better
  // than 'epsilon'.
  // NOTE: The number of trial steps is limited by 'fMaxNoSteps'. Integration will 
  //       stop if this maximum is reached, and the return value will be false.
  // On return
  //  - 'yOutput' provides the values at the end of the integration interval;
  //  - the return value is 'true' if integration succeeded to the end of the interval,
  //    and 'false' otherwise.

  #define PARTDEBUG

  typedef typename vecgeom::kVc::precision_v Double_v;
  typedef typename vecgeom::kVc::bool_v      Bool_v;
  // typedef vecgeom::Vector3D<Double_v>       ThreeVector_v;

#ifdef PARTDEBUG
  bool verbose = true;  
  if (partDebug) {
    if( verbose ) std::cout << " AccurateAdvance called with hStepIn= " ;
    for (int i = 0; i < nTracks; ++i) { std::cout<< hStepIn[i]<< " "; }
    std::cout<< std::endl;
  }
#endif 
  constexpr int ncompSVEC = TemplateGUFieldTrack<vecgeom::kVc>::ncompSVEC; // 12? - from TemplateGUFieldTrack
  // double numSteps[nTracks] = -1;

  // 'Invariants' during integration -- i.e. change only for a new track 
  Double_v  startCurveLength, chargeLane;
  Double_v  hTotalLane;      //  Total requested integration length, i.e.  x2 = x1 + hTotalLane
  Double_v  x1,  x2;         //  Start, end of integration interval (per lane)
  
  // Working values for lanes - can change at each iteration
  Double_v  y   [ncompSVEC], 
            dydx[ncompSVEC];
  Double_v  x;                   // Current value 
  Double_v  h,   hnext;          // Current & proposed-next step size
  Bool_v    succeededLane(true);
  Bool_v    isDoneLane(false);   // set true when there is a return statement 

  // Used together with KeepStepping only
  Double_v hDone(0.); // To keep track of length done in KeepStepping

  // Bool_v   lastStepSucceeded;  //  Not needed
  std::fill_n( succeeded, nTracks, 1); 

#ifdef COLLECT_STATISTICS
  // G4ThreadLocal  
  static int  noGoodSteps =0 ;  // Bad = chord > curve-len 
#endif 

  SetNTracks(nTracks);
  InitializeLanes( yInput, hStepIn, chargeIn, y, hTotalLane, chargeLane, startCurveLength );

  int trackNextInput = 4; 

  //  Ensure that hstep > 0
  if(!vecgeom::IsEmpty(hTotalLane<=0))
  {
    for (int i = 0; i < kVectorSize; ++i)
    {
      if (hTotalLane[i]<=0.0)
      {
        if (hTotalLane[i]==0.0)
        {
          std::cerr << "Proposed step is zero; hstep = " << hTotalLane[i] << " !";
          //succeededLane already true
        }
        else
        {
          std::cerr << "Invalid run condition." << std::endl
                  << "Proposed step is negative; hstep = " << hTotalLane[i] << std::endl;
          succeededLane[i] = false;
          succeeded    [i] = false; //the final bool array to be returned 
        }

        isDoneLane[i] = true;

        if (trackNextInput < nTracks)
        {
          // Insert new track because no processing required for h<=0 case 
          InsertNewTrack(yInput, hStepIn,      chargeIn,    i, trackNextInput, succeeded,
                         y,      hTotalLane, chargeLane, startCurveLength );
        }
      }
    }
  }
  
  x1= startCurveLength;
  x2= x1 + hTotalLane;   // x2 is also lane specific

  h = hTotalLane;        //  First integration length -- would be better to choose it carefully

  x = x1;

  // Recall strting value of y = ystart
  // for (i=0; i<fNoVars; i++)  { ystart[i] = y[i]; }

  Bool_v   lastStep(false);
  Double_v nstp(1); 

  // Double_v StartPosAr[3];

  // isDoneLane needed. In end, other conditions might keep changing
  // even if processing for that lane is finished. Need a way to store
  // the fact that the lane is done.
  // Either make a new single condition that combines isDoneLane 
  // and all other conditions or some conditions at least
  // For now, just adding isDoneLane : needs to be && or || with the first 3 
  // and keep nTracks condition in final ||
  // Say continue if isDoneLane is not 1111 and rest all conditions are not 0000
  // while ( !vecgeom::IsEmpty((nstp<=fMaxNoSteps) && (x < x2) && (!lastStep)) || trackNextInput < nTracks  )
  while ( (!vecgeom::IsFull(isDoneLane) && 
           !vecgeom::IsEmpty((nstp<=fMaxNoSteps) && (x < x2) && (!lastStep)) ) || 
           trackNextInput < nTracks  )
  {
#ifdef DEBUG
    std::cout<<"----hTotalLane is: "<< hTotalLane << std::endl;
#endif 
    // StartPosAr[0] = y[0]; StartPosAr[1] = y[1];  StartPosAr[2] = y[2]; 

    // Perform the Integration
    // -  Note:  Sequential code has "if( h > fMinimumStep )  {  ..." .  Concern ?

    Double_v hdid;
    Double_v numStepsLane; // Stat -- later can be used to avoid lane going on too long
    
    if (oneStep)
    {
      fpStepper->RightHandSideVIS( y, chargeLane, dydx );   // TODO: change to inline
      OneStep(      y,    dydx, chargeLane, epsilon, x, h, hdid, hnext);
      // Bool_v   lastStepSucceeded = (hdid == h);   
    }
    else
    {
       KeepStepping( y,    dydx, chargeLane, x, h, hTotalLane, epsilon, // lastLanes,
                     hnext, hDone, numStepsLane) ;
      // Differences from OneStep:
      //    hTotalLane   - 
      //    hDone        -
      //    numStepsLane - for 'score-keeping'
    }
    fNoTotalSteps++;

    // Ananya: discuss. What exactly is happening here?
    // h<=0 case: first condition false, second condition always true assuming smallest fraction and 
    // startCurveLength are positive. But what if startCurveLength is 0? Ask John what would happen
    // here for h<=0. : Ananya
    // If below bool is always true for h<=0 --> lastStep is true, hence the lane will be sent to 
    // StoreOutput.
    Bool_v avoidTinySteps = (h < epsilon * hTotalLane) || (h < fSmallestFraction * startCurveLength);
    lastStep = avoidTinySteps || lastStep;

    Double_v diff2 = (x2 - x);
    // For rest, check the proposed next stepsize 
    h = vecgeom::Max(hnext, fMinimumStep);
  #ifdef PARTDEBUG
    if (partDebug)
    {
      std::cout<< "diff before both masked assign statements is: "<<diff2 << std::endl;
      std::cout<< "h after checking proposed next stepsize is max. of : "<< hnext <<" and "<< x2-x << std::endl; 
      // std::cout.precision(17);  
      std::cout<<" Here, x2 is: "<<x2 << " and x is : "<< x << std::endl;
      Bool_v diffXAndX2 = (x+h) > x2 ;
      std::cout<< "Bool_v x+h > x2 is: "<< diffXAndX2 << std::endl;
    }
     Double_v diff = (x2 - x);
  #endif 

    // Ensure that the next step does not overshoot
    // vecgeom::MaskedAssign( x+h > x2, x2 - x, &h);
    h = vecgeom::Min( x2 -x , h);

  #ifdef PARTDEBUG
    if (partDebug)
    {
      Double_v hforDebug = x2 - x;
      std::cout<< "x2 -x is :  "<< hforDebug << std::endl;
      std::cout<< "diff is: " << diff << std::endl;
    }
  #endif 
  #ifdef PARTDEBUG
    if (partDebug)
    {
      // h = x2 - x; 
      std::cout<<"AccurateAdvance: hnext is: "<<hnext<<" and h is : "<<h<<std::endl; 
    }
  #endif 
    // When stepsize overshoots, decrease it!
    // Must cope with difficult rounding-error issues if hstep << x2

    lastStep = (h==0) || lastStep ;
#ifdef DEBUG
    if (partDebug)
    {
      std::cout<<" lastStep : "<< lastStep << std::endl;
    }
#endif
    nstp++;

    Bool_v CondNoOfSteps     = nstp<=fMaxNoSteps;
    Bool_v CondXLessThanx2   = x < x2;
    Bool_v CondIsNotLastStep = !lastStep; // lastStep is false

    bool condNoOfSteps     = vecgeom::IsFull(CondNoOfSteps    );
    bool condXLessThanx2   = vecgeom::IsFull(CondXLessThanx2  ); 
    bool condIsNotLastStep = vecgeom::IsFull(CondIsNotLastStep);

    Bool_v finishedLane;

    succeededLane = (x>=x2); // If it was a "forced" last step

    if ( !( condNoOfSteps && condXLessThanx2 && condIsNotLastStep ) )
    // Condition inside if can be stored in a variable and used for while condition. 
    // Saves some evaluations
    {
      finishedLane =  
              ( !CondNoOfSteps || !CondXLessThanx2 || !CondIsNotLastStep );
 #ifdef DEBUG
      if (partDebug)
      {
        std::cout<<" finishedLane:     "<< finishedLane     << std::endl;
        std::cout<<" CondNoOfSteps:    "<< CondNoOfSteps    << std::endl;
        std::cout<<" CondXLessThanx2:  "<< CondXLessThanx2  << std::endl;
        std::cout<<" CondIsNotLastStep:"<< CondIsNotLastStep<< std::endl;        
      }
#endif
      for (int i = 0; i < kVectorSize; ++i)
      {
        if (finishedLane[i] == 1 &&  fIndex[i] != -1)
        {
          // can be replaced with succeeded[fIndex[i]] = x[i] >= x2[i], one Vc vector reduced thus
          succeeded[fIndex[i]] = succeededLane[i]; //Final succeeded bool //  absorb into StoreOutput ??
          
          // Keep StoreOutput next allows succeeded to be changed (needed in case of h<=0 )
          // - We do not want to pass 'succeeded' to StoreOutput (why?) ... so don't absorb
          StoreOutput( y, x, yOutput, i, hStepIn, succeeded );

          // numSteps[ fIndex[i] ] = numStepsLane[i];
          
          if (trackNextInput<nTracks)
          {
            isDoneLane[i] = InsertNewTrack( yInput, hStepIn, chargeIn, i, trackNextInput,
                                            succeeded, y, hTotalLane, chargeLane, startCurveLength ); 

            // Logically the following should be part of InsertNewTrack().  
            //    Kept separate to reduce the number of its arguments.
            nstp     [i] = 1;
            lastStep [i] = false;
            x        [i] = x1[i];               // ?? Needed? Find something to set x < x2 
            h        [i] = hTotalLane[i];       // Needed ? Is it different ?
            x2       [i] = x[i] + hTotalLane[i];

            hDone[i] = 0.;         // It's new!
            numStepsLane[i]= 0;
          }
          else
          {
            isDoneLane[i] = true;
            fIndex    [i] = -1;
          }
        }
      }
    }
#ifdef DEBUG
    if (partDebug)
    {
      std::cout<<"Value of lastStep is: "<< lastStep <<std::endl;
      std::cout<<"isDoneLane is:        "<< isDoneLane <<std::endl;
    }
#endif 
       
/*    Bool_v leftLanes = (nstp<=fMaxNoSteps) && (x < x2) && (!lastStep) ;
    int countLeftLanes=0;
    int indLastLane;
    // std::cout << " leftLanes is: " << leftLanes << std::endl;
    if( !vecgeom::IsEmpty(leftLanes) )
    {
      for (int i = 0; i < kVectorSize; ++i)
      {
        if (leftLanes[i] == 1)
        {
          countLeftLanes++;
          indLastLane = i;
          // std::cout << indLastLane << std::endl;
        }
      }
    }

    // std::cout<< "countLeftLanes is: "<<countLeftLanes << std::endl;

    if (countLeftLanes == 1)
    {
      // double hstepOneLane = hTotalLane[indLastLane] - hDone[indLastLane];
      vecgeom::Vector3D<double> Pos, Mom;
      for (int i = 0; i < 3; ++i)
       {
         Pos[i] = y[i][indLastLane];
         Mom[i] = y[i+3][indLastLane];
       } 
      GUFieldTrack y_input(Pos, Mom); 
      GUFieldTrack y_output(Pos, Mom);
      // y_input.SetCurveLength( hDone[indLastLane] ) ;
      fpScalarDriver->AccurateAdvance(y_input, hStepIn[ fIndex[indLastLane] ] - hDone[indLastLane], epsilon, y_output );

      isDoneLane[indLastLane] == true;
      // Store Output
      double y_output_arr[12];
      y_output.DumpToArray(y_output_arr);
      yOutput[fIndex[indLastLane]].LoadFromArray(y_output_arr);
    }*/

  } // end of while loop

}  // end of AccurateAdvance ...........................
#endif /*NEWACCURATEADVANCE*/


// New constructor for KeepStepping method 
// Scalar stepper passed 
template <class Backend>
TemplateGUIntegrationDriver<Backend>
  ::TemplateGUIntegrationDriver( double  hminimum, 
                                 TemplateGUVIntegrationStepper<Backend> *pStepper,
                                 GUVIntegrationStepper                  *pScalarStepper,
                                 int     numComponents,
                                 int     statisticsVerbose                             )
  : TemplateGUIntegrationDriver( hminimum, 
                                 pStepper,
                                 numComponents,
                                 statisticsVerbose)
{
  fpScalarStepper = pScalarStepper; 
}

// New constructor. Takes in a scalar driver as well
template <class Backend>
TemplateGUIntegrationDriver<Backend>
  ::TemplateGUIntegrationDriver( double  hminimum, 
                                 TemplateGUVIntegrationStepper<Backend> *pStepper,
                                 GUVIntegrationStepper                  *pScalarStepper,
                                 GUIntegrationDriver                     *pScalarDriver,
                                 int     numComponents,
                                 int     statisticsVerbose                             )
  : TemplateGUIntegrationDriver( hminimum, 
                                 pStepper,
                                 numComponents,
                                 statisticsVerbose)
{
  fpScalarStepper = pScalarStepper; 
  fpScalarDriver  = pScalarDriver;
}


template <class Backend>
void
TemplateGUIntegrationDriver<Backend>
  ::SetPartDebug(bool debugValue)
{
  partDebug = debugValue;
}

template <class Backend>
void
TemplateGUIntegrationDriver<Backend>
  ::SetSteppingMethod(bool steppingMethod)
{
  oneStep = steppingMethod;
}

template <class Backend>
void
TemplateGUIntegrationDriver<Backend>
  ::ReportStepUnderflow( const std::string method, double x, double xnew, double hOriginal,
                          double hLast, double hNext,
                          double errPosSq, double errMomSq, const double dydx[], int numComponents)
{
   double errPos = ( errPosSq > 0.0 ) ? std::sqrt( errPosSq ) : errPosSq;
   double errMom = ( errMomSq > 0.0 ) ? std::sqrt( errMomSq ) : errMomSq;    
   std::cerr << "*** Warning from IntegrationDriver method " << method << std::endl
             << "    Stepsize *underflow* occurred " << std::endl ;
   std::cerr << "       Step's start x= " << x << " and end x= " << xnew << " coincide. "
             << std::endl
             << "       Last step size          = " << hLast << std::endl      
             << "       Proposed next step size = " << hNext << std::endl
             << "       Original (input) step   = " << hOriginal  << std::endl
             << "       Errors:  position = " << errPos
             <<                " momentum = " << errMom << std::endl
             << "       dy/dx : " << std::endl;
   for( int i=0; i < numComponents; ++i )
      std::cerr << "         [ " << i << " ] : " << dydx[i] << std::endl;
}

#endif /* TemplateGUIntegrationDriver_Def */
