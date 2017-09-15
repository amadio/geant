//
// Driver for explicit Runge-Kutta methods
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and ensures
//  that the errors are within acceptable bounds.
// When multiple tracks are integrated, provide different ways to 
//   handle the early end of some 'lanes' - while others continue.
//
// History:
//
// Adaptations of template interface: J. Apostolakis,     Nov 2017
// First templated version:  Ananya, Feb/March 2016
//     ( commit 95e1316bcc156a04c876d6ea0fc9e60a15eeac4f )
// 
// Adapted from G4MagInt_Drv class of Geant4 (G4MagIntegratorDriver)
//
// - Contributors: Ananya, J.Apostolakis                    2015-2017
// --------------------------------------------------------------------

#ifndef FlexibleIntegrationDriver_Def
#define FlexibleIntegrationDriver_Def

#include "TemplateFieldTrack.h"
#include "base/AlignedBase.h"
#include "FieldTrack.h"

// #include "TemplateGUVIntegrationStepper.h"
// #include "IntegrationStepper.h"

#include "base/Vector.h"

// Adding because adding scalar stepper for new constructor (KeepStepping)
#include "GUVIntegrationStepper.h"

// Adding to send in scalar driver to deal with 1/2 remaining lanes
#include "GUIntegrationDriver.h"
// #include "ScalarFieldTrack.h"

// #define NEWACCURATEADVANCE

// Eventual config: 
// template <class T_Stepper, unsigned int Nvar>

template <class Backend, class T_Stepper, unsigned int Nvar>
   class FlexibleIntegrationDriver : public vecgeom::AlignedBase
{

  public:

    template <typename T>
       using Vector3D = vecgeom::Vector3D<T>;

    using Real_v          = Backend;             // First adaptation
    // using Double_v        = Geant::Double_v;
    using Bool_v          = vecCore::Mask_v<Real_v>;
    using ThreeVectorSimd = Vector3D<Real_v>;
  
    FlexibleIntegrationDriver( double     hminimum,  //same 
                               T_Stepper *pStepper,
                               int        numberOfComponents  = 6,
                               int        statsVerbosity = 1     );

    FlexibleIntegrationDriver( const FlexibleIntegrationDriver& );
    // Copy constructor used to create Clone method
    ~FlexibleIntegrationDriver();

    // template <class Backend>  // Eventual config
    Bool_v  AccurateAdvance( const TemplateFieldTrack<Backend>& yInput,
                                   Real_v    hstep,
                                   Real_v    eps,                  // Requested y_err/hstep
                                   TemplateFieldTrack<Backend>& yOutput   );                            
                                   // Real_v  hinitial=0.0   );  // Suggested 1st interval
                                                                   // (e.g. from previous step) 
       // Above drivers for integrator (Runge-Kutta) with stepsize control. 
       // Integrates ODE starting values yInput 
       // from current s (s=s0) to s=s0+h with accuracy eps. 
       // On output yOutput is value at end of interval.
       // The concept is similar to the odeint routine from NRC 2nd edition p.721

     void  OneStepEach(   Real_v  ystart, // Like old RKF45step()
                    const Real_v  dydx,
                          Real_v& x,
                          Real_v  htry,
                          Real_v  eps,
                          Real_v& hdid,
                          Real_v& hnext ) ;
    
#ifdef NEWACCURATEADVANCE
     // 2nd AccurateAdvance
     // involves track insertion etc 
     // succeeded[] of length nTracks
     void  AccurateAdvance( /*const*/ FieldTrack yInput[],
                                  double      hstep[],
                                  double      epsilon,
                                  FieldTrack  yOutput[],
                                  int         nTracks, 
                                  bool        succeeded[] );

     void  OneStep(       Real_v  ystart[], // Like old RKF45step()
                    const Real_v  dydx[],
                          Real_v& x,
                          Real_v  htry,
                          Real_v  eps,      //  memb variables ?
                          Real_v& hdid,
                          Real_v& hnext ) ;

     void  InitializeAccurateAdvance( /*const*/ FieldTrack yInput[],
                                     const double      hstep [],
                                           Real_v    y[],
                                           Real_v    &hStepLane,
                                           Real_v    &startCurveLength);

     // Returns isDoneThisLane value i.e. whether the lane is done or not
     // True if the last lane left is an h<=0 lane, false for all else
     bool  InsertNewTrack( /*const*/ FieldTrack yInput[],
                          const double     hstep[],
                          const int        currIndex,
                                int&       trackNextInput,
                                bool       succeeded[],
                                Real_v   y[],
                                Real_v   &hStepLane,
                                Real_v   &startCurveLength );

     void  StoreOutput( const Real_v   y[],
                        const Real_v   x,
                              FieldTrack yOutput[],
                              int        currIndex,
                              double     hstep[],
                              bool       succeeded[]);

     void  SetNTracks( int nTracks );

     void  KeepStepping(       Real_v  ystart[], // Like old RKF45step()
                               Real_v  dydx[],
                               Real_v& x,
                               Real_v  htry,
                               Real_v  eps,      //  memb variables ?
                               Real_v& hdid,
                               Real_v& hnext,
                         const Real_v  hStepLane,
                               Real_v& hTotalDoneSoFar ) ;

     FlexibleIntegrationDriver( double hminimum,  //same 
                                T_Stepper *pStepper,
                                GUVIntegrationStepper   *pScalarStepper,
                                int    numberOfComponents  = 6,
                                int    statsVerbosity      = 1 );

     FlexibleIntegrationDriver( double hminimum,  //same 
                                T_Stepper             *pStepper,
                                GUVIntegrationStepper *pScalarStepper,
                                GUIntegrationDriver   *pScalarDriver,
                                int    numComponents  = 6,
                                int    statsVerbosity = 1      );

     void SetPartDebug (bool debugValue);

     // true for OneStep, false for KeepStepping
     void SetSteppingMethod (bool steppingMethod);
#endif

     Bool_v  QuickAdvance( TemplateFieldTrack<Backend>& y_posvel,        // INOUT
                           const Real_v      dydx[],  
                                 Real_v      hstep,           // IN
                           // Real_v&     dchord_step, //take out
                                 Real_v&     dyerr_pos_sq,
                                 Real_v&     dyerr_mom_rel_sq ) ;
       // New QuickAdvance that also just tries one Step
       //    (so also does not ensure accuracy)
       //    but does return the errors in  position and
       //        momentum (normalised: Delta_Integration(p^2)/(p^2) )

     // FlexibleIntegrationDriver* Clone() const;
       // Create an independent copy of the current object -- including independent 'owned' objects
       // NOTE: Evaluate whether this method is needed - 2017.11.09
     
     // VEquationOfMotion<Backend>* GetEquationOfMotion() { return fpStepper->GetEquationOfMotion(); }
     // const VEquationOfMotion<Backend>* GetEquationOfMotion() const { return fpStepper->GetEquationOfMotion(); } 
     
     // Auxiliary methods
     inline double GetHmin()        const { return fMinimumStep; } 
     inline double GetSafety()      const { return fSafetyFactor; }
     inline double GetPowerShrink() const { return fPowerShrink; }
     inline double GetPowerGrow()   const { return fPowerGrow; } 
     inline double GetErrcon()      const { return fErrcon; }
     
     inline void   GetDerivatives( const TemplateFieldTrack<Backend> &y_curr,     // const, INput
                                        Real_v    charge, 
                                        Real_v    dydx[]   );  //       OUTput
        // Accessors.

     // ===> Set methods suppressed - parameters now const 2017.11.08
     //
     
     // inline void RenewStepperAndAdjust(T_Stepper *Stepper);
        // Sets a new stepper 'Stepper' for this driver. Then it calls
        // ReSetParameters to reset its parameters accordingly.

     // inline void ReSetParameters(double new_safety= 0.9 ); //same 
        //  i) sets the exponents (fPowerGrow & fPowerShrink), 
        //     using the current Stepper's order, 
        // ii) sets the safety
        // ii) calculates "fErrcon" according to the above values.

     // inline void SetSafety(double valS) ;       // {fSafetyFactor= valS; }
     // inline void SetPowerShrink(double valPs) ; // { fPowerShrink= valPs; } 
     // inline void SetPowerGrow (double valPg) ;  // { fPowerGrow= valPg; }
     // inline void SetErrcon(double valEc) ;      // { fErrcon= valEc; }
        // When setting safety or fPowerGrow, errcon will be set to a 
        // compatible value.

     inline double ComputeAndSetErrcon();

     inline const T_Stepper* GetStepper() const;
     inline       T_Stepper* GetStepper();

     void  OneGoodStep(       Real_v  ystart[], // Like old RKF45step()
                        const Real_v  dydx[],
                              Real_v& x,
                              Real_v  htry,
                              Real_v  eps,      //  memb variables ?
                              Real_v& hdid,
                              Real_v& hnext ) ;
        // This takes one Step that is as large as possible while 
        // satisfying the accuracy criterion of:
        // yerr < eps * |y_end-y_start|

     Real_v ComputeNewStepSize( Real_v  errMaxNorm,    // normalised error
                                  Real_v  hstepCurrent); // current step size
        // Taking the last step's normalised error, calculate
        // a step size for the next step.
        // Do not limit the next step's size within a factor of the
        // current one.

     Real_v ComputeNewStepSize_WithinLimits( Real_v  errMaxNorm,    // normalised error
                                               Real_v  hstepCurrent); // current step size
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
     void WarnSmallStepSize( Real_v , // hnext, 
                             Real_v , // hstep, 
                             Real_v , // h,     
                             Real_v , // xDone,
                             int    /*noSteps*/) {} //;  //warnings per track, probably neeed to change all to double :ananya
     /* void WarnTooManySteps ( Real_v , x1start, 
                             Real_v x2end, 
                             Real_v xCurrent) {} //; */
     /* void WarnEndPointTooFar( Real_v  endPointDist, 
                              Real_v  hStepSize , 
                              Real_v  epsilonRelative,
                              int     debugFlag) {}  */
     //  Issue warnings for undesirable situations

  /***********
     //add index in order to print one at a time :ananya 
     void PrintStatus(  const Real_v*      StartArr,
                              Real_v       xstart,
                        const Real_v*      CurrentArr, 
                              Real_v       xcurrent, 
                              Real_v       requestStep, 
                              int          subStepNo ) {} //;
     void PrintStatus(  const TemplateFieldTrack<Backend>&  StartFT,
                        const TemplateFieldTrack<Backend>&  CurrentFT, 
                              double       requestStep, 
                              int          subStepNo ) {} //;
     void PrintStat_Aux( const TemplateFieldTrack<Backend>& aScalarFieldTrack,
                               double      requestStep, 
                               double      actualStep,
                               int         subStepNo,
                               double      subStepSize,
                               double      dotVelocities ) {} //;       
  ***********/

     //  Verbose output for debugging

     void PrintStatisticsReport() {} //;
     //  Report on the number of steps, maximum errors etc.

#ifdef QUICK_ADV_ARRAY_IN_AND_OUT      
     Bool_v QuickAdvance(      Real_v     yarrin[],     // In
                         const Real_v     dydx[],  
                               Real_v     hstep,        
                               Real_v     yarrout[],    // Out
                               Real_v&    dchord_step,  // Out
                               Real_v&    dyerr );      // in length
#endif

   private:

     FlexibleIntegrationDriver& operator=(const FlexibleIntegrationDriver&);
     // Private copy constructor and assignment operator.

   private:

     // ---------------------------------------------------------------
     // DEPENDENT Objects
     T_Stepper *fpStepper;
     
     // ---------------------------------------------------------------
     //  INVARIANTS 

     double  fMinimumStep; //same 
     // Minimum Step allowed in a Step (in absolute units)
     double  fSmallestFraction; //same      //   Expected range 1e-12 to 5e-15;  
     // Smallest fraction of (existing) curve length - in relative units
     //  below this fraction the current step will be the last 

     const int  fNoIntegrationVariables;  // Number of Variables in integration
     const int  fMinNoVars;               // Minimum number for TemplateFieldTrack<Backend>
     const int  fNoVars;                  // Full number of variable

     int   fMaxNoSteps;
     static const int  fMaxStepBase;  

     const double fSafetyFactor;
     const double fPowerShrink;   //  exponent for shrinking
     const double fPowerGrow;    //  exponent for growth
     const double fErrcon;
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
     Real_v fDyerr_max, fDyerr_mx2;
     Real_v fDyerrPos_smTot, fDyerrPos_lgTot, fDyerrVel_lgTot; 
     Real_v fSumH_sm, fSumH_lg; 
     // Step Statistics 

     int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
     // Could be varied during tracking - to help identify issues
#ifdef NEWACCURATEADVANCE
     //Variables required for track insertion algorithm
     constexpr unsigned int kVectorSize = 4; // TODO: template it on the backend 
     int fIndex[kVectorSize]; // int *fIndex; 
     int fNTracks;
     int fStepperCalls = 0;
     GUVIntegrationStepper *fpScalarStepper;
     GUIntegrationDriver *fpScalarDriver;
     bool partDebug = false ; 
     bool oneStep = true; // false for KeepStepping
#endif 
};

// #include "GUIntegratorDriver.icc"

template <class Backend, class T_Stepper, unsigned int Nvar>
constexpr double FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>::fMaxSteppingIncrease;

template <class Backend, class T_Stepper, unsigned int Nvar>
constexpr double FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>::fMaxSteppingDecrease;


template <class Backend, class T_Stepper, unsigned int Nvar>
inline
double FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
   ::ComputeAndSetErrcon()
{
  fErrcon = std::pow(fMaxSteppingIncrease/fSafetyFactor,1.0/fPowerGrow);
  return fErrcon;
} 

// START Methods to modify parameters ----- Suppressed
#ifdef NON_CONST_PARAMS
template <class Backend, class T_Stepper, unsigned int Nvar>
inline
void FlexibleIntegrationDriver<Backend, T_Stepper, unsigned int Nvar>
  ::ReSetParameters(double new_safety)
{
  fSafetyFactor = new_safety;
  fPowerShrink  = -1.0 / fpStepper->IntegratorOrder();
  fPowerGrow    = -1.0 / (1.0 + fpStepper->IntegratorOrder());
  ComputeAndSetErrcon();
}

template <class Backend, class T_Stepper, unsigned int Nvar>
inline
void FlexibleIntegrationDriver<Backend, T_Stepper, unsigned int Nvar>
  ::SetSafety(double val)
{ 
  fSafetyFactor=val;
  ComputeAndSetErrcon();
}

template <class Backend, class T_Stepper, unsigned int Nvar>
inline
void FlexibleIntegrationDriver<Backend, T_Stepper, unsigned int Nvar>
  ::SetPowerGrow(double  val)
{ 
  fPowerGrow=val;
  ComputeAndSetErrcon(); 
}

template <class Backend, class T_Stepper, unsigned int Nvar>
inline
void FlexibleIntegrationDriver<Backend, T_Stepper, unsigned int Nvar>
  ::SetErrcon(double val)
{ 
  fErrcon=val;
}

template <class Backend, class T_Stepper, unsigned int Nvar>
inline
void FlexibleIntegrationDriver<Backend, T_Stepper, unsigned int Nvar>
  ::RenewStepperAndAdjust(T_Stepper *pStepper)
{  
  fpStepper = pStepper; 
  ReSetParameters();
}
#endif
// END:   Methods to modify parameters ----- Suppressed

template <class Backend, class T_Stepper, unsigned int Nvar>
inline
const T_Stepper* FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::GetStepper() const
{
  return fpStepper;
}

template <class Backend, class T_Stepper, unsigned int Nvar>
inline
T_Stepper* FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::GetStepper() 
{
  return fpStepper;
}

template <class Backend, class T_Stepper, unsigned int Nvar>
inline
int FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::GetMaxNoSteps() const
{
  return fMaxNoSteps;
}

template <class Backend, class T_Stepper, unsigned int Nvar>
inline
void FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::SetMaxNoSteps(int val)
{
  fMaxNoSteps= val;
}

template <class Backend, class T_Stepper, unsigned int Nvar>
inline
void FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::GetDerivatives(const TemplateFieldTrack<Backend> &y_curr, // const, INput
                         Real_v  charge, 
                         Real_v  dydx[])  // OUTput
{ 
  Real_v  tmpValArr[TemplateFieldTrack<Backend>::ncompSVEC];
  y_curr.DumpToArray( tmpValArr  );
  fpStepper -> RightHandSideVIS( tmpValArr , charge, dydx );
}


template <class Backend, class T_Stepper, unsigned int Nvar>
const int  FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>::fMaxStepBase = 250;  // Was 5000

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
template <class Backend, class T_Stepper, unsigned int Nvar>
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::FlexibleIntegrationDriver( double  hminimum, 
                                 T_Stepper *pStepper,
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
  // fIndex = new int[kVectorSize];
}


//  Copy Constructor - used by Clone
template <class Backend, class T_Stepper, unsigned int Nvar>
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
::FlexibleIntegrationDriver( const FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>& right ) 
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
   const T_Stepper *protStepper = right.GetStepper();
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
template <class Backend, class T_Stepper, unsigned int Nvar>
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::~FlexibleIntegrationDriver()
{ 
  if( fStatisticsVerboseLevel > 1 )
  {
    PrintStatisticsReport();
  }

  // delete[] fIndex;
  // delete fpScalarDriver;
  // delete fpScalarStepper;
  // delete fpStepper;
}

#if 0
template <class Backend, class T_Stepper, unsigned int Nvar>
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>* 
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::Clone() const
{
   return new FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>(*this);
}
#endif

// ---------------------------------------------------------
template <class Backend, class T_Stepper, unsigned int Nvar>
void
   FlexibleIntegrationDriver<double,T_Stepper,Nvar>::OneGoodStep(  double y[],        // InOut
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

  int ncompSVEC = TemplateFieldTrack<double>::ncompSVEC;
  double yerr[ncompSVEC], ytemp[ncompSVEC];

  bool verbose= false;
  if( verbose ) 
     std::cout << "OneGoodStep called with htry= " << htry << std::endl;
  
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
    fpStepper-> StepWithErrorEstimate(y,dydx,h,ytemp,yerr);
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

// ----------------------------------------------------------------------------------

template <class Backend, class T_Stepper, unsigned int Nvar>
bool 
FlexibleIntegrationDriver<double>  // scalar
  ::AccurateAdvance(const TemplateFieldTrack<double>& yInput,  // scalar
                          double                      hstep,
                          double                      epsilon,
                          TemplateFieldTrack<double>& yOutput )
                          // Real_v  hinitial)
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

  typedef vecgeom::Vector3D<double>   ThreeVector;

  // std::cout<<"\n " << std::endl;
#ifdef COLLECT_STATISTICS
  constexpr double perMillion  = 1.0e-6;
  constexpr double perThousand = 1.0e-3;
  int no_warnings = 0;
#endif 

  int nstp, i;
  double x, hnext, hdid, h;

  constexpr int ncompSVEC = TemplateFieldTrack<double>::ncompSVEC; // 12 ? 

#ifdef GUDEBUG_FIELD
  static int dbg=1;
  static int nStpPr=50;   // For debug printing of long integrations
  double ySubStepStart[ncompSVEC];

  // std::cout << " AccurateAdvance called with hstep= " << hstep << std::endl;
#endif

  double  y     [ncompSVEC], 
            dydx  [ncompSVEC];
  double  ystart[ncompSVEC], 
            yEnd  [ncompSVEC]; 
  double  x1, x2;
  Bool_v succeeded(true), lastStepSucceeded;

  double startCurveLength;

  int  noFullIntegr=0, noSmallIntegr = 0 ;
  // G4ThreadLocal
  static int  noGoodSteps =0 ;  // Bad = chord > curve-len 
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


  double StartPosAr[3];
  double charge(-1.);

  while ( ((nstp++)<=fMaxNoSteps) && (x < x2) && (!lastStep) )
  {
    StartPosAr[0]= y[0]; StartPosAr[1]= y[1];  StartPosAr[2]= y[2]; 
    
    fpStepper->RightHandSideVIS( y, charge, dydx );   // TODO: change to inline

#ifdef COLLECT_STATISTICS
    fNoTotalSteps++;
#endif 

    // Perform the Integration
    if( h > fMinimumStep )
    { 
      OneGoodStep(y,dydx,x,h,epsilon,hdid,hnext) ;

      //--------------------------------------
      lastStepSucceeded= (hdid == h);   

    }

#ifdef COLLECT_STATISTICS
    if (lastStepSucceeded)  { noFullIntegr++; }
    else                    { noSmallIntegr++; }
#endif 

    ThreeVector EndPos( y[0], y[1], y[2] );

    // Check the endpoint
    const double edx= y[0] - StartPosAr[0];
    const double edy= y[1] - StartPosAr[1];
    const double edz= y[2] - StartPosAr[2];
    double endPointDist2= vecgeom::Sqrt(edx*edx+edy*edy+edz*edz) ; // (EndPos-StartPos).Mag(); 

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

/*template <class Backend, class T_Stepper, unsigned int Nvar>
void
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
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
  std::cerr << " WARNING from FlexibleIntegrationDriver::WarnSmallStepSize():  " ; // << std::endl;
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
  // G4Exception("FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>::WarnSmallStepSize()", "GeomField1001",
  //             JustWarning, message);
  noWarningsIssued++;
}

// ---------------------------------------------------------

template <class Backend, class T_Stepper, unsigned int Nvar>
void
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::WarnTooManySteps( double x1start, 
                      double x2end, 
                      double xCurrent)
{

   std::cerr << "WARNING from FlexibleIntegrationDriver::WarnTooManySteps()" << std::endl;
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
   // G4Exception("FlexibleIntegrationDriver::WarnTooManySteps()", "GeomField1001",
   //             JustWarning, message);
}

// ---------------------------------------------------------

template <class Backend, class T_Stepper, unsigned int Nvar>
void
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
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
    std::cerr << "WARNING in FlexibleIntegrationDriver::WarnEndPointTooFar()" << std::endl;
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
    // G4Exception("FlexibleIntegrationDriver::WarnEndPointTooFar()", "GeomField1001",
    //             JustWarning, message);
  }
}*/

// ---------------------------------------------------------

template </*class Backend*/>
void
FlexibleIntegrationDriver<<Real_v, T_Stepper, Nvar>>
  ::OneGoodStep(       Real_v y[],        // InOut
                 const Real_v dydx[],
                       Real_v& x,         // InOut
                       Real_v htry,
                       Real_v eps_rel_max,
                       Real_v& hdid,      // Out
                       Real_v& hnext      )    // Out

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
  Real_v errmax_sq;
  Real_v h, htemp, xnew ;

  Real_v yerr [TemplateFieldTrack<<Real_v, T_Stepper, Nvar>>::ncompSVEC], 
           ytemp[TemplateFieldTrack<<Real_v, T_Stepper, Nvar>>::ncompSVEC];

  // std::cout << "OneGoodStep called with htry= " << htry << std::endl;
  
  h = htry ; // Set stepsize to the initial trial value

  Real_v inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  Real_v errpos_sq=0.0;    // square of displacement error
  Real_v errmom_sq=0.0;    // square of momentum vector difference

  int iter;

  static int tot_no_trials=0;  // thread_local
  const int max_trials=100; 

  for (iter=0; iter<max_trials ;iter++)
  {
    if (true)
    {
      tot_no_trials++;
      fpStepper-> StepWithErrorEstimate(y,dydx,h,ytemp,yerr);

      Real_v eps_pos = eps_rel_max * vecgeom::Max(h, fMinimumStep);  // Uses remaining step 'h'
      Real_v inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos);

      // Evaluate accuracy
      errpos_sq =  yerr[0]*yerr[0] + yerr[1]*yerr[1] + yerr[2]*yerr[2] ; 
      errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

      // Accuracy for momentum
      Real_v magmom_sq=  y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
      Real_v sumerr_sq=  yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5]; 

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
      Real_v errPower = Vc::exp( (0.5*fPowerShrink)*vecgeom::Log(errmax_sq) );
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
  Real_v errPower = Vc::exp( (0.5*GetPowerGrow())* vecgeom::Log(errmax_sq));
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
template <class Backend, class T_Stepper, unsigned int Nvar>//
typename Mask<Backend>
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::QuickAdvance( TemplateFieldTrack<Backend>&       y_posvel,         // INOUT
                  const Real_v  dydx[],  
                        Real_v  hstep,       // In
                  // Real_v& dchord_step,
                        Real_v& dyerr_pos_sq,
                        Real_v& dyerr_mom_rel_sq )  
{
  // typedef typename Backend Real_v;
//  typedef typename Mask<Backend>      Bool_v;
  // Real_v dyerr_pos_sq, dyerr_mom_rel_sq;  
  Real_v yerr_vec[TemplateFieldTrack<Backend>::ncompSVEC],
           yarrin  [TemplateFieldTrack<Backend>::ncompSVEC], 
           yarrout [TemplateFieldTrack<Backend>::ncompSVEC]; 
  Real_v s_start;
  Real_v dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

  static int no_call=0;  // thread_local 
  no_call ++; 

  // Move data into array
  y_posvel.DumpToArray( yarrin );      //  yarrin  <== y_posvel 
  s_start = y_posvel.GetCurveLength();

  // Do an Integration Step
  fpStepper-> StepWithErrorEstimate(yarrin, dydx, hstep, yarrout, yerr_vec) ; 
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
template <class Backend, class T_Stepper, unsigned int Nvar>
typename Mask<Backend>  
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::QuickAdvance(       Real_v     yarrin[],    // In
                  const Real_v     dydx[],  
                        Real_v     hstep,       // In
                        Real_v     yarrout[],
                        Real_v&    dchord_step,
                        Real_v&    dyerr )      // In length
{
   std::cerr << "ERROR in FlexibleIntegrationDriver::QuickAdvance()" << std::endl;
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
template <class Backend, class T_Stepper, unsigned int Nvar>
Real_v 
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::ComputeNewStepSize( Real_v  errMaxNorm,    // max error  (normalised)
                        Real_v  hstepCurrent)  // current step size
{
  Real_v hnew;

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
template <class Backend, class T_Stepper, unsigned int Nvar>
Real_v 
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::ComputeNewStepSize_WithinLimits( Real_v  errMaxNorm,    // max error  (normalised)
                                     Real_v  hstepCurrent)  // current step size
{
  Real_v hnew, htemp;

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
/*template <class Backend, class T_Stepper, unsigned int Nvar>
void 
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
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
   ScalarFieldTrack  StartFT(ThreeVector(0,0,0),
                 ThreeVector(0,0,0), 0., 0. );
   ScalarFieldTrack  CurrentFT (StartFT);

   StartFT.LoadFromArray( StartArr, fNoIntegrationVariables); 
   StartFT.SetCurveLength( xstart);
   CurrentFT.LoadFromArray( CurrentArr, fNoIntegrationVariables); 
   CurrentFT.SetCurveLength( xcurrent );

   PrintStatus(StartFT, CurrentFT, requestStep, subStepNo ); 
}

// ---------------------------------------------------------------------------
template <class Backend, class T_Stepper, unsigned int Nvar>
void FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::PrintStatus( const TemplateFieldTrack<Backend>&  StartFT,
                 const TemplateFieldTrack<Backend>&  CurrentFT, 
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
              << " FlexibleIntegrationDriver: Current Position  and  Direction" << " "
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
template <class Backend, class T_Stepper, unsigned int Nvar>
void FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::PrintStat_Aux( const TemplateFieldTrack<Backend>&  aScalarFieldTrack,
                   double             requestStep, 
                   double             step_len,
                   int                subStepNo,
                   double             subStepSize,
                   double             dotVeloc_StartCurr)
{
    const ThreeVector Position=      aScalarFieldTrack.GetPosition();
    const ThreeVector UnitVelocity=  aScalarFieldTrack.GetMomentumDirection();
 
    if( subStepNo >= 0)
    {
       std::cout << std::setw( 5) << subStepNo << " ";
    }
    else
    {
       std::cout << std::setw( 5) << "Start" << " ";
    }
    double curveLen= aScalarFieldTrack.GetCurveLength();
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
    // std::cout << std::setw( 7) << aScalarFieldTrack.GetKineticEnergy();
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
template <class Backend, class T_Stepper, unsigned int Nvar>
void 
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::PrintStatisticsReport()
{
  int noPrecBig= 6;
  int oldPrec= std::cout.precision(noPrecBig);

  std::cout << "FlexibleIntegrationDriver Statistics of steps undertaken. " << std::endl;
  std::cout << "FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>: Number of Steps: "
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
}*/
 
// ---------------------------------------------------------------------------
template <class Backend, class T_Stepper, unsigned int Nvar>
void FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
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
template <class Backend, class T_Stepper, unsigned int Nvar>
void 
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::SetNTracks(int nTracks)
// Set fNTracks 
{
  fNTracks = nTracks;
}

// #define DEBUG

template </*class Backend*/>
void 
FlexibleIntegrationDriver<<Real_v, T_Stepper, Nvar>>
  ::InitializeAccurateAdvance(/*const*/ FieldTrack yInput[],
                              const double     hstep [],
                                    Real_vy[],
                                    Real_v&hStepLane,
                                    Real_v&startCurveLength)
// Initialization step for AccurateAdvance/
// Converts input scalar stream to acceptable form of Vc vectors
// for vector processing in OneStep
{
#ifdef DEBUG
  std::cout<<"----Initializing AccurateAdvance----"<<std::endl;
#endif

  double yStartScalar[fNoVars]; //fNoVars or fNoIntegrationVariables? Ask : Ananya
  for (int j = 0; j < kVectorSize; ++j)
  {
    fIndex          [j] = j;
    hStepLane       [j] = hstep[j];
    yInput          [j].DumpToArray(yStartScalar);
    startCurveLength[j] = yInput[j].GetCurveLength();
    for (int i = 0; i < fNoVars; ++i)
    {
      y[i][j] = yStartScalar[i]; 
    }
  }

} // End of InitializeAccurateAdvance function


template </*class Backend*/>
bool 
FlexibleIntegrationDriver<<Real_v, T_Stepper, Nvar>>
  ::InsertNewTrack( /*const*/ FieldTrack                    yInput[],
                    const double                        hstep[],
                    const int                           currIndex,
                          int&                          trackNextInput,
                          bool                          succeeded[],
                          Real_vy[],
                          Real_v&hStepLane,
                          Real_v&startCurveLength )
// Inserts a new track whenever a lane is finished. 
// returns isDoneLane = true for h<=0 case, false otherwise
// because in former case, no further work is required
{
#ifdef DEBUG
  std::cout<<"----Inserting New Track "<< trackNextInput << " at position "<< currIndex <<std::endl;
#endif 

  bool done = true; // to get the while loop starting
  while(trackNextInput < fNTracks && done== true)
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
      done = false;
      double yScalar[fNoVars];
      yInput[trackNextInput].DumpToArray(yScalar);
      for (int i = 0; i < fNoVars; ++i)
      {
        y[i][currIndex] = yScalar[i];
      }
      fIndex          [currIndex] = trackNextInput;
      hStepLane       [currIndex] = hstep[trackNextInput];
      startCurveLength[currIndex] = yInput[trackNextInput].GetCurveLength();

    }

    trackNextInput++;
    
  }
  
  return done;

} // End of InsertNewTrack function


template </*class Backend*/>
void 
FlexibleIntegrationDriver<<Real_v, T_Stepper, Nvar>>
  ::StoreOutput(const Real_vy[],
                const Real_vx,
                      FieldTrack   yOutput[],
                      int          currIndex,
                      double       hstep[],
                      bool         succeeded[] )
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




template </*class Backend*/>
void
FlexibleIntegrationDriver<<Real_v, T_Stepper, Nvar>>
  ::OneStep(       Real_v y[],        // InOut
             const Real_v dydx[],
                   Real_v& x,         // InOut
                   Real_v htry,
                   Real_v eps_rel_max,
                   Real_v& hdid,      // Out
                   Real_v& hnext      )    // Out
// Derived from OneGoodStep
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
  if (partDebug)
  {
    std::cout<<"\n"<<std::endl;  
  }
#endif 

  Real_v errmax_sq;
  Real_v h, htemp, xnew ;

  int ncompSVEC = TemplateFieldTrack<Real_v>::ncompSVEC;

  Real_v yerr [ncompSVEC], 
           ytemp[ncompSVEC];

  h = htry ; // Set stepsize to the initial trial value

  Real_v inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  Real_v errpos_sq=0.0;    // square of displacement error
  Real_v errmom_sq=0.0;    // square of momentum vector difference

  int iter;

  static int tot_no_trials=0;  // thread_local
  const int max_trials=100; 

  int finished[kVectorSize] = {0,0,0,0};

  Real_v hFinal, hnextFinal, xFinal, hdidFinal, errmax_sqFinal;
  Real_v yFinal[ncompSVEC];
  Bool_v errMaxLessThanOne(false), hIsZeroCond(false);

  for (iter = 0; iter < max_trials; iter++)
  {
    if (  !vecgeom::IsFull(hIsZeroCond || errMaxLessThanOne) )
    {
      tot_no_trials++;
      fpStepper-> StepWithErrorEstimate(y,dydx,h,ytemp,yerr);
      fStepperCalls++;

      if (0)
      {
        std::cout<< "----h is: " << h[0] << " at iter: " << iter << std::endl;
        // std::cout<< " yerr is: " << yerr[0] << std::endl;
      }

#ifdef PARTDEBUG
      if (0)
      {
        std::cout<< "----yerr is: " << yerr[0][0] <<" "<<yerr[1][0]<<" "<<yerr[2][0] << " " << yerr[3][0] << " " << yerr[4][0] << " " << yerr[5][0]  << std::endl;
      }
#endif 
      Real_v eps_pos = eps_rel_max * vecgeom::Max(h, fMinimumStep);  // Uses remaining step 'h'
      Real_v inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos);

      // Evaluate accuracy
      errpos_sq =  yerr[0]*yerr[0] + yerr[1]*yerr[1] + yerr[2]*yerr[2] ; 
      errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

      // Accuracy for momentum
      Real_v magmom_sq=  y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
      Real_v sumerr_sq=  yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5]; 

      vecgeom::CondAssign(magmom_sq > 0.0, sumerr_sq/magmom_sq, sumerr_sq, &errmom_sq);


      errmom_sq *= inv_eps_vel_sq;
      errmax_sq = vecgeom::Max( errpos_sq, errmom_sq ); // Square of maximum error
#ifdef PARTDEBUG
      if (0)
      {
        std::cout<< "----eps_pos is: "<< eps_pos[0] << std::endl;
        std::cout<< "----inv_eps_pos_sq is: "<< inv_eps_pos_sq[0] << std::endl;
        std::cout<< "----errmom_sq is: "<< errmom_sq[0] << std::endl;
        std::cout<< "----errpos_sq is: "<< errpos_sq[0] << std::endl;
        std::cout<< "----errmax_sq is: "<< errmax_sq[0] << std::endl;
      }

#endif 
      errMaxLessThanOne = ( errmax_sq <=1.0 );
      if ( !vecgeom::IsEmpty(errMaxLessThanOne) )
      {
        for (int i = 0; i < kVectorSize; ++i)
        {
          // Probably could use several MaskedAssigns as well
          if ( errMaxLessThanOne[i] ==1 && finished[i] !=-1 )
          {
            /* StoreFinalValues() */
            finished      [i] = -1;
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

      // Step failed; compute the size of retrial Step.
      // Ananya : adding a statement. Later check the sanity or work around
      Real_v errPower = Vc::exp( (0.5*fPowerShrink)*vecgeom::Log(errmax_sq) ); 
      htemp = fSafetyFactor *h* errPower;
      // htemp = fSafetyFactor *h* vecgeom::Pow( errmax_sq, 0.5*fPowerShrink );
      // Can use the loop below instead of the lines above since power is 
      // expensive operation. 
/*      for (int i = 0; i < kVectorSize; ++i)
      {
        if (finished[i] != -1)
        {
          htemp[i] = fSafetyFactor *h[i]* std::pow(errmax_sq[i], 0.5*fPowerShrink);
        }
      }*/

      h = vecgeom::Max(htemp, 0.1*h);
      
      xnew = x + h;

      hIsZeroCond = (xnew == x);
      if ( !vecgeom::IsEmpty(hIsZeroCond) )
      {
        for (int i = 0; i < kVectorSize; ++i)
        {
          // Probably could use several MaskedAssigns as well
          if ( hIsZeroCond[i] ==1 && finished[i] !=-1 )
          {
            /* StoreFinalValues() */
            finished      [i] = -1;
            hFinal        [i] = h[i];
            errmax_sqFinal[i] = errmax_sq[i];
            for (int j = 0; j < ncompSVEC; ++j)
            {
              yFinal[j][i] = ytemp[j][i];
            }
          }
        }
      }      
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
#ifdef PARTDEBUG
  if (partDebug)
  {
    std::cout << "FlexibleIntDrv: 1--step - Loop done at iter = " << iter << " with htry= " << htry <<std::endl;
  }
#endif 

  h         = hFinal;
  errmax_sq = errmax_sqFinal;

  // Compute size of next Step
  Real_v errPower = Vc::exp( (0.5*GetPowerGrow())* vecgeom::Log(errmax_sq));
  hnext = GetSafety() * h * errPower;
  // hnext = GetSafety()*vecgeom::Pow(errmax_sq, 0.5*GetPowerGrow());
  vecgeom::MaskedAssign(errmax_sq <= fErrcon*fErrcon, fMaxSteppingIncrease*h, &hnext); // No more than a factor of 5 increase

  x += (hdid = h);

  for(int k=0;k<fNoIntegrationVariables;k++) { y[k] = yFinal[k]; }

#ifdef PARTDEBUG
  if (partDebug)
  {
    std::cout<< " hdid= "<<hdid<<" and hnext= "<<hnext<<  std::endl;
  }
#endif
  return;
}   // end of  OneStep .............................


template </*class Backend*/>
void
FlexibleIntegrationDriver<<Real_v, T_Stepper, Nvar>>
  ::KeepStepping(       Real_v  y[],        // InOut
                        Real_v  dydx[],
                        Real_v& x,         // InOut
                        Real_v  htry,
                        Real_v  eps_rel_max,
                        Real_v& hdid,      // Out
                        Real_v& hnext,     // Out
                  const Real_v  hStepLane,
                        Real_v& hTotalDoneSoFar  )    

// Derived from OneGoodStep
// WIP

{
#ifdef PARTDEBUG
  if (partDebug)
  {
    std::cout<<"\n"<<std::endl;
  }
#endif 

  Real_v errmax_sq;
  Real_v h, htemp, xnew ;

  Real_v yerr [TemplateFieldTrack<Real_v>::ncompSVEC], 
           ytemp[TemplateFieldTrack<Real_v>::ncompSVEC];

  h = htry ; // Set stepsize to the initial trial value

  Real_v inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

  Real_v errpos_sq = 0.0;    // square of displacement error
  Real_v errmom_sq = 0.0;    // square of momentum vector difference

  int iter;


  static int tot_no_trials = 0;  // thread_local
  const  int max_trials    = 100; 

  int finished[kVectorSize] = {0}; // This makes all elements of array 0

  Real_v hFinal(0.), hnextFinal, xFinal, hdidFinal, errmax_sqFinal;
  Real_v yFinal[TemplateFieldTrack<Real_v>::ncompSVEC]; // = y[]
  for (int i = 0; i < TemplateFieldTrack<Real_v>::ncompSVEC; ++i)
  {
    yFinal[i] = y[i];
  }

  Bool_v errMaxLessThanOne(false), hIsZeroCond(false);

  Real_v x2 = x + (hStepLane - hTotalDoneSoFar);
  // Real_v x2 = x + htry;
  Bool_v  errMaxLessThanOneLocal(false), hIsZeroCondLocal(false);
  // int htryExhausted[kVectorSize] = {0};
  Bool_v   htryExhausted(false);
  Real_v charge(+1.);

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
      fpStepper-> StepWithErrorEstimate( yFinal, dydx, h, ytemp, yerr);
      fStepperCalls++;

#ifdef DEBUG
      if (partDebug)
      {
        std::cout<< "\n----yerr is: " << yerr[0] <<" "<<yerr[1]<<" "<<yerr[2] << std::endl;
      }
#endif 
      Real_v eps_pos = eps_rel_max * vecgeom::Max(h, fMinimumStep);  // Uses remaining step 'h'
      Real_v inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos);

      // Evaluate accuracy
      errpos_sq  =  yerr[0]*yerr[0] + yerr[1]*yerr[1] + yerr[2]*yerr[2] ; 
      errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

      // Accuracy for momentum
      Real_v magmom_sq =  y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
      Real_v sumerr_sq =  yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5]; 

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
            finished         [i] = -1;
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
              for (int j = 0; j < TemplateFieldTrack<Real_v>::ncompSVEC; ++j)
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
              hTotalDoneSoFar[i] += h[i]; // Diff. from hFinal because can be carried over 
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
      Real_v errPower = Vc::exp( (0.5*fPowerShrink)*vecgeom::Log(errmax_sq) ); 
      htemp = GetSafety() *h* errPower;


      h = vecgeom::Max(htemp, 0.1*h);

      h = vecgeom::Min(h, hStepLane - hTotalDoneSoFar);
      
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
            finished   [i] = -1;
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
              for (int j = 0; j < TemplateFieldTrack<Real_v>::ncompSVEC; ++j)
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
              hTotalDoneSoFar[i] += h[i]; // Diff. from hFinal because can be carried over 
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
              for (int j = 0; j < TemplateFieldTrack<Backend>::ncompSVEC; ++j)
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
  Real_v errPower = Vc::exp( (0.5*GetPowerGrow())* vecgeom::Log(errmax_sq));
  hnext = GetSafety() * h * errPower;
  vecgeom::MaskedAssign(errmax_sq <= fErrcon*fErrcon, fMaxSteppingIncrease*h, &hnext); // No more than a factor of 5 increase

  // std::cout<<"fPowerShrink is: "<<1/fPowerShrink<<" and fPowerGrow is: "<<1/GetPowerGrow()<< std::endl;

  x += (hdid = h);

  for(int k=0;k<fNoIntegrationVariables;k++) { y[k] = yFinal[k]; }

#ifdef PARTDEBUG
  if (partDebug)
  {
    std::cout << "FlexibleIxntDrv: 1--step - Loop done at iter = " << iter << " with htry= " << htry <<std::endl;
    std::cout<< " hdid= "<<hdid<<" and hnext= "<<hnext<<  std::endl;
    std::cout<< "htryExhausted is: "<< htryExhausted << std::endl;
  }
#endif 

  return;
}   // end of  KeepStepping .............................


#ifdef NEWACCURATEADVANCE
template </*class Backend*/>
void
FlexibleIntegrationDriver<<Real_v, T_Stepper, Nvar>>
  ::AccurateAdvance(/*const*/ FieldTrack yInput[],
                          double     hstep[],
                          double     epsilon,
                          FieldTrack yOutput[],
                          int        nTracks,
                          bool       succeeded[])
{
  // Built on original AccurateAdvance. Takes buffer stream of nTracks
  // Converts them to Vc vectors for processing
  // Inserts new track when processing for a lane is finished.

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

  #define PARTDEBUG

  using ThreeVector = vecgeom::Vector3D<Real_v>;

  Real_v x, hnext, hdid, h;

  int ncompSVEC = TemplateFieldTrack<Real_v>::ncompSVEC; //12, to be derived from TemplateFieldTrack

#ifdef GUDEBUG_FIELD
  // static int dbg=1;
  // static int nStpPr=50;   // For debug printing of long integrations
  // Real_v ySubStepStart[ncompSVEC];
#endif

#ifdef PARTDEBUG
  if (partDebug)
  {
    std::cout << " AccurateAdvance called with hstep= " ;
    for (int i = 0; i < nTracks; ++i)
    {
      std::cout<< hstep[i]<< " ";
    }
    std::cout<< std::endl;
  }
#endif 

  Real_v  y   [ncompSVEC], 
            dydx[ncompSVEC];
  Real_v  x1, x2;
  std::fill_n( succeeded, nTracks, 1); 

  Bool_v lastStepSucceeded;

  Real_v startCurveLength;

  // G4ThreadLocal
#ifdef COLLECT_STATISTICS
  static int  noGoodSteps =0 ;  // Bad = chord > curve-len 
#endif 

  Real_v hStepLane; 
  Real_v hTotalDoneSoFar(0.); // To keep track of hDone in KeepStepping
  Bool_v   succeededLane(true); 
  Bool_v   isDoneLane(false); // set true when there is a return statement 
  int trackNextInput = 4; 

  SetNTracks(nTracks);
  InitializeAccurateAdvance( yInput, hstep, y, hStepLane, startCurveLength );

  //  Ensure that hstep > 0
  if(!vecgeom::IsEmpty(hStepLane<=0))
  {
    for (int i = 0; i < kVectorSize; ++i)
    {
      if (hStepLane[i]<=0.0)
      {
        if (hStepLane[i]==0.0)
        {
          std::cerr << "Proposed step is zero; hstep = " << hstep << " !";
          //succeededLane already true
        }
        else
        {
          std::cerr << "Invalid run condition." << std::endl
                  << "Proposed step is negative; hstep = " << hstep << "." << std::endl;
          succeededLane[i] = false;
          succeeded    [i] = false; //the final bool array to be returned 
        }

        isDoneLane[i] = true;

        if (trackNextInput < nTracks)
        {
    // Insert new track because no processing required for h<=0 case 
          InsertNewTrack(yInput, hstep, i, trackNextInput, succeeded, y, hStepLane, startCurveLength );
        }
      }
    }
  }

  x1= startCurveLength; 
  x2= x1 + hStepLane; // x2 also needs to be lane specific

  h = hStepLane; // = x2 -x1 ; or x2 - x

  x = x1;

  // Why both y and ystart? Ask John : Ananya
  // for (i=0; i<fNoVars; i++)  { y[i] = ystart[i]; }

  Bool_v   lastStep(false);
  Real_v nstp(1); 

  Real_v StartPosAr[3];

  // Ananya : making random charge now
  // needs to be passed in some other way finally
  Real_v charge(-1.);

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
    std::cout<<"----hStepLane is: "<< hStepLane << std::endl;
#endif 
    StartPosAr[0] = y[0]; StartPosAr[1] = y[1];  StartPosAr[2] = y[2]; 
    

    // Perform the Integration
    // Ananya: doing OneGoodStep always. Ask John what should be done.
    // if( h > fMinimumStep )
    // { 

    if (oneStep)
    {
      fpStepper->RightHandSideVIS( y, charge, dydx );   // TODO: change to inline
      OneStep( y, dydx, x, h, epsilon, hdid, hnext);
    }
    else KeepStepping( y, dydx, x, h, epsilon, hdid, hnext, hStepLane, hTotalDoneSoFar) ;

    fNoTotalSteps++;
    // KeepStepping( y, dydx, x, h, epsilon, hdid, hnext, hStepLane, hTotalDoneSoFar) ;
    lastStepSucceeded = (hdid == h);   
    // }

    ThreeVector EndPos( y[0], y[1], y[2] );

    // Check the endpoint
    const Real_v edx= y[0] - StartPosAr[0];
    const Real_v edy= y[1] - StartPosAr[1];
    const Real_v edz= y[2] - StartPosAr[2];
    Real_v endPointDist2= vecgeom::Sqrt(edx*edx+edy*edy+edz*edz) ; 

    // Ananya: discuss. What exactly is happening here?
    // h<=0 case: first condition false, second condition always true assuming smallest fraction and 
    // startCurveLength are positive. But what if startCurveLength is 0? Ask John what would happen
    // here for h<=0. : Ananya
    // If below bool is always true for h<=0 --> lastStep is true, hence the lane will be sent to 
    // StoreOutput.
    Bool_v avoidNumerousSmallSteps = (h < epsilon * hStepLane) || (h < fSmallestFraction * startCurveLength);
    lastStep = avoidNumerousSmallSteps || lastStep;

    Real_v diff2 = (x2 - x);
    // For rest, check the proposed next stepsize 
    h = vecgeom::Max(hnext, fMinimumStep);
  #ifdef PARTDEBUG
    if (0)
    {
      std::cout<< "diff before both masked assign statements is: "<<diff2 << std::endl;
      std::cout<< "h after checking proposed next stepsize is max. of : "<< hnext <<" and "<< x2-x << std::endl; 
      // std::cout.precision(17);  
      std::cout<<" Here, x2 is: "<<x2 << " and x is : "<< x << std::endl;
      Bool_v diffXAndX2 = (x+h) > x2 ;
      std::cout<< "Bool_v x+h > x2 is: "<< diffXAndX2 << std::endl;
    }
     Real_v diff = (x2 - x);
  #endif 


    // Ensure that the next step does not overshoot
    // vecgeom::MaskedAssign( x+h > x2, x2 - x, &h);
    h = vecgeom::Min(x2 -x , h);

  #ifdef PARTDEBUG
    if (0)
    {
      Real_v    hforDebug = x2 - x;
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
          succeeded[fIndex[i]] = succeededLane[i]; //Final succeeded bool // might be absorbed in StoreOutput
          
          // Keep StoreOutput after succeeded[fIndex[i]] = succeededLane[i]; so that
          // succeeded can be changed again. Needs to be changed in case of h<=0
          // If succeeded is completely absorbed in StoreOutput, then succeededLane also needs
          // to be passed, which we do not want to do. 
          StoreOutput( y, x, yOutput, i, hstep, succeeded );

          if (trackNextInput<nTracks)
          {
            isDoneLane[i] = InsertNewTrack( yInput, hstep, i, trackNextInput, 
                                            succeeded, y, hStepLane, startCurveLength ); 

            nstp      [i] = 1;     // logically part of InsertNewTrack, not done so to reduce
            lastStep  [i] = false; // number of parameters to be passed to the function
            x         [i] = x1[i]; // ?? Needed? Find something to set x<x2 
            h         [i] = hStepLane[i];// Can absorb in InsertNewTrack as well, leads to too many variables though
                                        // Maybe ask John interpretation of this h and then put in InsertNewTrack
                                        // with appropriate name 
            x2        [i] = x[i] + hStepLane[i];

            hTotalDoneSoFar[i] = 0.; // Setting to 0 for every new track inserted. 
                                     // Adding here so as not to pollute InsertNewTrack
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
      // double hstepOneLane = hStepLane[indLastLane] - hTotalDoneSoFar[indLastLane];
      vecgeom::Vector3D<double> Pos, Mom;
      for (int i = 0; i < 3; ++i)
       {
         Pos[i] = y[i][indLastLane];
         Mom[i] = y[i+3][indLastLane];
       } 
      ScalarFieldTrack y_input(Pos, Mom); 
      ScalarFieldTrack y_output(Pos, Mom);
      // y_input.SetCurveLength( hTotalDoneSoFar[indLastLane] ) ;
      fpScalarDriver->AccurateAdvance(y_input, hstep[ fIndex[indLastLane] ] - hTotalDoneSoFar[indLastLane], epsilon, y_output );

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
template <class Backend, class T_Stepper, unsigned int Nvar>
FlexibleIntegrationDrivertemplate <Backend, T_Stepper, Nvar>   
  ::FlexibleIntegrationDriver( double                  hminimum, 
                               T_Stepper              *pStepper,
                               GUVIntegrationStepper  *pScalarStepper,
                               int                     numComponents,
                               int                     statsVerbose )
  : FlexibleIntegrationDriver( hminimum, 
                                 pStepper,
                                 numComponents,
                               statsVerbose),
   fpScalarStepper(pScalarStepper)
{
}

// New constructor. Takes in a scalar driver as well
template <class Backend, class T_Stepper, unsigned int Nvar>
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::FlexibleIntegrationDriver( double  hminimum, 
                                 T_Stepper *pStepper,
                                 GUVIntegrationStepper                  *pScalarStepper,
                                 GUIntegrationDriver                     *pScalarDriver,
                                 int     numComponents,
                                 int     statisticsVerbose                             )
  : FlexibleIntegrationDriver( hminimum, 
                                 pStepper,
                                 numComponents,
                                 statisticsVerbose)
{
  fpScalarStepper = pScalarStepper; 
  fpScalarDriver  = pScalarDriver;
}


template <class Backend, class T_Stepper, unsigned int Nvar>
void
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::SetPartDebug(bool debugValue)
{
  partDebug = debugValue;
}

template <class Backend, class T_Stepper, unsigned int Nvar>
void
   
FlexibleIntegrationDriver<Backend, T_Stepper, Nvar>
  ::SetSteppingMethod(bool steppingMethod)
{
  oneStep = steppingMethod;
}

#endif /* FlexibleIntegrationDriver_Def */
