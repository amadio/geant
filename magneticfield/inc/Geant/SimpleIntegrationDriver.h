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

#ifndef SimpleIntegrationDriver_Def
#define SimpleIntegrationDriver_Def

// #include "Geant/TemplateFieldTrack.h"
#include "base/AlignedBase.h"
#include "Geant/FieldTrack.h"

#include "base/Vector.h"

#include "Geant/FlexIntegrationDriver.h"

#include "Geant/BaseRkIntegrationDriver.h"
#include "Geant/AuxVecMethods.h"

#include "ErrorEstimatorSixVec.h"
    
#include "Geant/FormattedReporter.h"  // Is direct include needed ?

#include "Geant/PrintDriverProgress.h"

#include "Geant/VectorTypes.h" //  Defines geant::Double_v
#include "Geant/math_wrappers.h"

#define CHECK_ONE_LANE 1
//  Define to check a single lane

#define CONST_DEBUG    1
//  Define to turn 'partDebug' into compile time constant

#define DRIVER_PRINT_PROGRESS   1

#ifdef CHECK_ONE_LANE
#include "IntegrationDriverConstants.h"
#endif

// --------------------------------------------------------------
template <class T_Stepper, unsigned int Nvar>
   class SimpleIntegrationDriver : public BaseRkIntegrationDriver<T_Stepper, Nvar>, 
                                       // FlexIntegrationDriver,
                                   public vecgeom::AlignedBase
{
public:
  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

  SimpleIntegrationDriver(double     hminimum, // same
                          T_Stepper *pStepper,
                          double     epsRelMax,
                          int        numberOfComponents = 6 );

  virtual ~SimpleIntegrationDriver();

  virtual void AccurateAdvance( const FieldTrack yInput[],
                                const double hstep[],
                                const double charge[], // double epsilon,
                                FieldTrack yOutput[],
                                int nTracks,
                                bool succeeded[]
                              ) const override final;

#ifdef EXTEND_SINGLE
  virtual void AccurateAdvance( const FieldTrack & yInput,
                                const double       hstep,
                                const double       charge, // double epsilon,
                                FieldTrack       & yOutput,
                                bool               succeeded
                              ) const override final;
#endif

  // Implemented in terms of the following templated function:
  template <class Real_v>
     void AccurateAdvance(const FieldTrack yInput[], const double hstep[], const double charge[], /* double epsilon, */
                       FieldTrack yOutput[], bool succeeded[], int nTracks) const;
  // Drive Runge-Kutta integration of ODE for several tracks (ntracks)
  // with starting values yInput, from current 's'=0 to s=h with variable
  // stepsize to control error, so that it is bounded by the relative
  // accuracy eps.  On output yOutput is value at end of interval.
  // The concept is similar to the odeint routine from NRC 2nd edition p.721

  // Auxiliary methods
  inline double GetHmin() const { return fMinimumStep; }
  // inline double GetSafetyFactor() const { return kSafetyFactor; }
  // inline constexpr double GetPowerShrink() const { return kPowerShrink; }
  // inline constexpr double GetPowerGrow() const { return kPowerGrow; }
  inline double GetErrcon() const { return fErrcon; }

  inline int GetMaxNoSteps() const { return fMaxNoSteps; }

  unsigned long GetNumberOfStepperCalls() { return fStepperCalls; }
  unsigned long GetNumberOfTotalSteps() { return fNoTotalSteps; } // fNumberOfTotalSteps; }
  /*****
      inline void   GetDerivatives( const TemplateFieldTrack<Real_v> &y_curr,     // const, INput
                                          Real_v    charge,
                                          Real_v    dydx[]   );  //       OUTput
   ******/

  // EquationOfMotion<Real_v>* GetEquationOfMotion() { return fpStepper->GetEquationOfMotion(); }
  // const EquationOfMotion<Real_v>* GetEquationOfMotion() const { return fpStepper->GetEquationOfMotion(); }

  // SimpleIntegrationDriver* Clone() const;
  // Create an independent copy of the current object -- including independent 'owned' objects
  // NOTE: Evaluate whether this method is needed - 2017.11.09

#ifdef QUICK_ADV_ARRAY_IN_AND_OUT
  template <class Real_v>
  vecCore::Mask_v<Real_v> QuickAdvance(Real_v yarrin[], // In
                                       const Real_v dydx[], Real_v hstep,
                                       Real_v yarrout[],    // Out
                                       Real_v &dchord_step, // Out
                                       Real_v &dyerr);      // in length
#endif

protected:
  // Implementation methods
  template <class Real_v>
  void OneGoodStep(const Real_v yStart[], //  [N]
                   const Real_v dydx[],   //  [N]
                   const Real_v charge,
                   Real_v &x, // InOut
                   Real_v htry,
                   double epsilon, // Was const Real_v  epsilon,  // required relative accuracy
                   Real_v yEnd[],  // [N]
                   Real_v &hdid, Real_v &hnext) const;
  // Integrate over the length htry 'together', until each lane succeeds
  //   if not for the initial lenght 'htry', then for a reduced size step
  //    ( reduced in order to likely succeed. )
  // In this version each lane stops as soon with its first success.

  template <class Real_v>
  int InitializeLanes(const FieldTrack yInput[], const double hstep[], const double charge[], int nTracks,
                      int      indexArr[], // [vecCore::VectorSize<Real_v>()] - Output
                      Real_v   y[],     // [Nvar]        - Output
                      Real_v & hStepLane, Real_v &chargeLane, Real_v &startCurveLength,
                      int    & numFilled, // [Out]: number of filled lanes
                      bool     badStepSize[],   // Output - needed to reset values
                      int    & numBadSize) const;
  // Load 'work' into array y, indices of lanes into indexArr, etc
  // Return array index of 'next' track - i.e. first track not yet loaded

  template <class Real_v>
  bool InsertNewTrack(const FieldTrack yInput[], const double hstep[], const double charge[], const int slot,
                      int    & trackNextInput, bool succeeded[], Real_v y[], Real_v &hStepLane, Real_v &chargeLane,
                      Real_v & startCurveLength,
                      int      indexArr[], // [vecCore::VectorSize<Real_v>()]
                      int      nTracks,
                      bool     badStepSize[], 
                      int    & numBadSize ) const;
     
  template <class Real_v>
  void StoreOutput(const Real_v yEnd[], const Real_v x,
                   int currIndex, // Index in Real_v
                   FieldTrack yOutput[],
                   int indOut, // location in Ouptut
                   const double hstep[], bool succeeded[], int nTracks) const;

  // void ComputeAndSetErrcon() { BaseRkIntegrationDriver<T_Stepper, Nvar>::ComputeAndSetErrcon(); }  // Compute dependent parameters
  // double GetErrcon() const {  return BaseRkIntegrationDriver<T_Stepper, Nvar>::GetErrcon(); }
  // double GetErrcon() const {  double val= BaseRkIntegrationDriver<T_Stepper,Nvar>::GetErrcon();
  //    std::cout << " <== SiD gets errcon = " << val << " ==> ";
  //    return val; 
  // }
  
  void ComputeAndSetErrcon(); 
  double GetErrcon() const { assert(fErrcon > 0.0); return  fErrcon; }
  
  void CheckParameters()     { BaseRkIntegrationDriver<T_Stepper, Nvar>::CheckParameters(); }
  int IncrementStepperCalls() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::IncrementStepperCalls() ; }
  int IncrementNumberSteps()  { return ++fNoTotalSteps; }
  // Setting parameters ( few now )

  // Compute dependent parameters

  // Check
  // void CheckParameters(); // Sanity check of

public: 
  // Access parameters
  // double GetPowerShrink()  const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fPowerShrink; }
  // double GetPowerGrow()    const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fPowerGrow; }
  double GetSafetyFactor() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fSafetyFactor; }
  int    GetVerboseLevel() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fVerboseLevel; }
  int    GetStatisticsVerboseLevel() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStatisticsVerboseLevel(); }
  double GetSmallestFraction() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetSmallestFraction(); }

  // Accessors.
  unsigned int GetStepperOrder() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStepperOrder(); }
  unsigned long GetNumberOfStepperCalls() { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetNumberOfStepperCalls(); }
  unsigned long GetNumberOfTotalSteps()   { return fNoTotalSteps; } 
     // BaseRkIntegrationDriver<T_Stepper,Nvar>::GetNumberOfTotalSteps(); }

  double GetMinimumStep() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetMinimumStep(); }

  inline unsigned int GetMaxNoSteps() const { return fMaxNoSteps; } 
  // inline unsigned int GetMaxNoSteps() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetMaxNoSteps() ; } // fMaxNoSteps; }  
  
  // Modifier methods
  void SetMaxNoSteps(unsigned int val) { fMaxNoSteps = val; }
  // void SetMaxNoSteps(int val) { BaseRkIntegrationDriver<T_Stepper,Nvar>::SetMaxNoSteps(val) ; }
  
  int IncrementStepperCalls()  { return BaseRkIntegrationDriver<T_Stepper,Nvar>::IncrementStepperCalls() ; }

  // const T_Stepper *GetStepper() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStepper(); }
  //       T_Stepper *GetStepper()       { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStepper(); }
  
  const T_Stepper *GetStepper() const { return fpStepper; }   // Tried suppressing to check its use  2019.04.10 @ 17:20
        T_Stepper *GetStepper() { return fpStepper; }

  void CreateInput(FieldTrack yInput[], int nTracks)
  { BaseRkIntegrationDriver<T_Stepper,Nvar>::CreateInput( yInput, nTracks);  }
  // void CreateInput(FieldTrack yInput[], int nTracks);
        
  // For debugging:
  
#ifdef CHECK_ONE_LANE
  mutable int    laneToCheck = -1; //  Lane that will be checked / printed, if any
#endif          
  
public: // For now
  template <class Real_v>
  bool TestInitializeLanes(); // (int numTracks)
                              // Simple, unit-test like check.  Returns ok or not
  template <class Real_v>
  bool CheckOutput(Real_v yOutput[], int lane, int initialSlot, std::string testName, std::string varName);
  // returns 'true' if ok

  void ReportInvalidStepInputs(double hStep[], int nTracks);

private:
  // Private copy constructor and assignment operator.

  SimpleIntegrationDriver(const SimpleIntegrationDriver &);
  // Copy constructor used to create Clone method

  SimpleIntegrationDriver &operator=(const SimpleIntegrationDriver &) = delete;

  template <class Real_v>
  void StoreGoodValues(const Real_v yWork[], const Real_v &hValue, const Real_v &epsSqVal,
                       const vecCore::Mask_v<Real_v> &storeFlag, Real_v yFinal[], Real_v &hFinalStore,
                       Real_v &epsSqStore) const;
  // Auxiliary method, to store results of selected ('good') lanes

private:
  // ---------------------------------------------------------------
  // DEPENDENT Objects
  T_Stepper *fpStepper;

  // ---------------------------------------------------------------
  //  INVARIANTS

  // Moved to BaseRkIntegrationDriver
  //
  // double fMinimumStep; // same
  //    Minimum Step allowed in a Step (in absolute units)
  // static constexpr double fSmallestFraction = 1.0e-7; // Expected value: larger than 1e-12 to 5e-15;
  //    Smallest fraction of (existing) curve length - in relative units
  //    below this fraction the current step will be the last

  const double fHalfPowerShrink;
  unsigned long fMaxNoSteps= 100;  // Default - expect it to be overriden in constructor
  
  // ---------------------------------------------------------------
  // Compilation constants
#ifdef CONST_DEBUG  
  const bool partDebug  = true;  // false;                 // Enforce debugging output
#endif
  const int ncompSVEC   = FieldTrack::NumCompFT; // expect 6, later 8, eventually up to 12
  const bool useOneStep = true;                  //  Algorithm selection - false for KeepStepping

  // ---------------------------------------------------------------
  //  INVARIANTS

  double fMinimumStep; // same
  // Minimum Step allowed in a Step (in absolute units)
  static constexpr double fSmallestFraction = 1.0e-7; // Expected value: larger than 1e-12 to 5e-15;
  // Smallest fraction of (existing) curve length - in relative units
  //  below this fraction the current step will be the last

  //  INVARIANTS - continued 

  // const int  fNoIntegrationVariables;  // Number of Variables in integration
  // const int fMinNoVars; // Minimum number for TemplateFieldTrack<Real_v>
  // const int fNoVars;    // Full number of variable

  static constexpr int fMaxStepBase = 250;

  // static constexpr double fSafetyFactor= 0.9; // -> Failed to compile on clang 9.1 2017.12.05
  static constexpr double kSafetyFactor = 0.9;                                            //     OK ...
  static constexpr double kPowerShrink  = -1.0 / T_Stepper::GetIntegratorOrder();         //  exponent for shrinking
  static constexpr double kPowerGrow    = -1.0 / (1.0 + T_Stepper::GetIntegratorOrder()); //  exponent for growth
  /*const*/ double fErrcon = 0.0;
  // Parameters used to grow and shrink trial stepsize.

  // double fSurfaceTolerance = 1.e-6;

  //  Stepsize can increase by no more than 5.0
  //           and decrease by no more than x10. = 0.1
  static constexpr double fMaxSteppingIncrease = 5.0;
  static constexpr double fMaxSteppingDecrease = 0.1;
  // Maximum stepsize increase/decrease factors.

  // int fStatisticsVerboseLevel         = 0;
  mutable unsigned long fStepperCalls = 0UL;
  // ---------------------------------------------------------------
  //  STATE
  // public:

  // Related to AccurateAdvance
  // Variables required for track insertion algorithm
  // int    fNTracks = 0;        //  number under current integration -> ~ current array size
  // ScalarIntegrationStepper *fpScalarStepper= nullptr;  // --> call Stepper with scalar values (args)
  // ScalarIntegrationDriver  *fpScalarDriver = nullptr;  // --> (later) use this driver with scalar args
#ifndef CONST_DEBUG    
  mutable bool partDebug = false ;
#endif  
  // ---

  mutable unsigned long fNoTotalSteps = 0, fNoBadSteps = 0, fNoSmallSteps = 0, fNoInitialSmallSteps = 0;

  int fVerboseLevel; // Verbosity level for printing (debug, ..)
                     // Could be varied during tracking - to help identify issues

}; // End of class definition -- SimpleIntegrationDriver

#ifdef DRIVER_PRINT_PROGRESS  
#include "Geant/FormattedReporter.h"
#endif

//  Constructor
//
template <class T_Stepper, unsigned int Nvar>
SimpleIntegrationDriver<T_Stepper, Nvar>::SimpleIntegrationDriver(double     hminimum,
                                                                  T_Stepper *pStepper,
                                                                  double     epsRelMax,                                                                  
                                                                  int        numComponents )
    :
      BaseRkIntegrationDriver<T_Stepper, Nvar>(hminimum,
                                               pStepper,     // or pass pStepper->GetIntegratorOrder()  ??
                                               epsRelMax,
                                               numComponents,
                                               false),  // statistics Verbosity
      fHalfPowerShrink( 0.5 * BaseRkIntegrationDriver<T_Stepper, Nvar>::GetPowerShrink() )
{
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
  assert(pStepper != nullptr);
  assert(Nvar <= (unsigned int)ncompSVEC); // Ensure that arrays are large enough for Integr.

  fpStepper = pStepper;

  // fMaxNoSteps = fMaxStepBase / fpStepper->GetIntegratorOrder();
  SetMaxNoSteps( fMaxStepBase / fpStepper->GetIntegratorOrder() );

  ComputeAndSetErrcon();
  std::cout << " Errcon = " << GetErrcon() << std::endl;
  assert( GetErrcon() > 0.0 );
  
  CheckParameters();

#ifdef GUDEBUG_FIELD
  fVerboseLevel = 2;
#endif

  if (fVerboseLevel) {
    std::cout << "SiD:ctor> Stepper Order= " << pStepper->GetIntegratorOrder() << " > Powers used: "
              << " shrink (half) = " << kPowerShrink     << "  grow (half) = " << kPowerGrow << std::endl;
              << " shrink (full) = " << GetPowerShrink() << "  grow (full) = " << GetPowerGrow() << std::endl
              << " and with  max-relative-error = " << epsRelMax << std::endl;     
  }
  if ((GetVerboseLevel() > 0) || (GetStatisticsVerboseLevel() > 1)) {
     std::cout << "SimpleIntegrationDriver created. " << std::endl;
     //         << "invE_nS, QuickAdv-2sqrt with Statistics " << fStatsStatus << std::endl;
     // ( fStatsEnabled ? "enabled" : "disabled" )
  }

  // For track insertion
}

// ---------------------------------------------------------

//  Copy Constructor - used by Clone
template <class T_Stepper, unsigned int Nvar>
SimpleIntegrationDriver<T_Stepper, Nvar>::SimpleIntegrationDriver(
    const SimpleIntegrationDriver<T_Stepper, Nvar> &right)
    :
   BaseRkIntegrationDriver<T_Stepper, Nvar>::BaseRkIntegrationDriver( right )
{
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
   
  // const T_Stepper *protStepper = right.GetStepper();
  // fpStepper                    = protStepper->Clone();
  // ==> This should happen in the base class .....
   
  ComputeAndSetErrcon();
  assert( GetErrcon() > 0.0 );
  
  SetMaxNoSteps( fMaxStepBase / GetStepper()->GetIntegratorOrder() );  

  if ((GetVerboseLevel() > 0) || (GetStatisticsVerboseLevel() > 1)) {
    std::cout << "SimpleIntegrationDriver copy Constructor. "
              << std::endl;
  }
}

// ---------------------------------------------------------

//  Destructor
template <class T_Stepper, unsigned int Nvar>
SimpleIntegrationDriver<T_Stepper, Nvar>::
   ~SimpleIntegrationDriver()
   //======================
{
   if (GetStatisticsVerboseLevel() > 1) {
      // PrintStatisticsReport();
   }
}

// ---------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void SimpleIntegrationDriver<T_Stepper, Nvar>::OneGoodStep(const Real_v   yStart[],
                                                           const Real_v   dydx[],
                                                           const Real_v   charge,
                                                           Real_v       & x,       // InOut  - ToDo: rename to xCurrent
                                                           Real_v         htry,
                                                           double         eps_rel_max,
                                                           // const Real_v  eps_rel_max,
                                                           Real_v         yFinal[], // Out-values
                                                           Real_v       & hdid,    // Out - achieved length
                                                           Real_v       & hnext)   // Out - proposed next integration length
   const
// This version:  J. Apostolakis,  13 November 2017.
//   Lanes are integrated until all have either,
//     - succeeded with the initial interval (1st iteration),
//     -     >>    at a later iteration, with a reduced step size, or
//     - failed due to step-size underlow.
//  That is, no minimum step size exists (or it is not respected.)
//
// A maximum of number of iterations is observed.
//
// Simplest method, meant as baseline, or for cases where all lanes are
//  expected to succeed in 1 step in a large fraction of cases.
// -------------------------------------------------------------------------
//
// Derived from OneGoodStep
//
// Driver for one Runge-Kutta Step with monitoring of local truncation error
// to ensure accuracy and adjust stepsize.
//  Input are dependent variable array y[0,...,5] and its derivative dydx[0,...,5]
// at the starting value of the independent variable x . Also input are stepsize
// to be attempted htry, and the required accuracy eps. On output y and x
// are replaced by their new values, hdid is the stepsize that was actually
// accomplished, and hnext is the estimated next stepsize.
// Similar to function rkqs from Numerical Recipes in C, 2nd Ed:p. 719
//

{
  using std::cout;  using std::cerr;  using std::endl;
  using Bool_v = vecCore::Mask_v<Real_v>;

  using vecCore::Get;
  using vecCore::math::Exp;
  using vecCore::math::Log;
  using vecCore::math::Max;
  using vecCore::math::Min;

  using FormattedReporter::ReportManyRowsOfDoubles;
  using FormattedReporter::ReportRowOfBools;
  using FormattedReporter::ReportRowOfDoubles;
  using FormattedReporter::ReportRowOfSquareRoots;
  using FormattedReporter::ReportRowOfDoublesIf;  
  using FormattedReporter::ReportOneLane;
  // using ReportValuesOfVectors::ReportConditionLanes;
  using PrintDriverProgress::ReportConditionLanes;
   
  // if (partDebug) { cout << "\n" << endl; }
  
  static  std::atomic<unsigned int> numCalls(0);
  int  currentCallNo= numCalls++;
  
  const Real_v xStart(x);
  Real_v xnew;
  Real_v yerr[Nvar], ytemp[Nvar];

  Real_v h = htry; // Set stepsize to the initial trial value
  // Rename it to hStep or hCurrent -- ie the current step size
#ifdef STORE_ONCE  
  Real_v errmaxSqFallThru(0.0);
#endif  
  static int tot_no_trials = 0; // Should be thread_local - or suppressed. Just statistics
  const int max_trials     = 100;

#ifdef CHECK_ONE_LANE
  const int trackToPrint = IntegrationDriverConstants::GetInstance()->GetTrackToCheck();
#ifndef CONST_DEBUG
  partDebug = ( laneToCheck != -1 );
#endif
  // bool  printNow = (laneToCheck >= 0) && vecCore::Get( htry, laneToCheck ) > 0.0 ;       
  if( laneToCheck != -1 ) {
     cout << "* SID:AccAdv Global Vars> trackToPrint = " << trackToPrint << " l2c / laneToCheck = " << laneToCheck;
     cout << "  Args>  h[l2c] = " << vecCore::Get( htry, laneToCheck );     
     cout << endl;
  }
#endif
  if( partDebug ) { std::cout << "SID: OneGoodStep called.   Call number " << currentCallNo << std::endl; }
  
  // int finishedArr[vecCore::VectorSize<Real_v>()] = {0,0,0,0};
  Bool_v finished = (htry <= 0.); //  Allows h <=0 as signal lane is empty. // Was = false;

  // vecCore::Int_v  finishedInt = 0;

  Real_v hFinal = htry, errmax_sqFinal = Real_v(0.); // xFinal, hdidFinal,
  // Real_v yFinal[ncompSVEC];
  Bool_v goodStep(false), stepSizeUnderflow(false);

  // for (iter = 0; iter < max_trials; iter++)
  unsigned int iter = 0;
  int itersLeft     = max_trials;

  // ReportManyRowsOfDoubles( "yStart",  yStart, Nvar );
  // Real_v magmomInit_sq = yStart[3] * yStart[3] + yStart[4] * yStart[4] + yStart[5] * yStart[5];

  do {
    Real_v errmax_sq = 0.0;

    // Bool_v alreadyFinished = finished;  // State at start of iteration
    Bool_v Active = !finished;

#ifdef CHECK_ONE_LANE    
    bool  printLane = (laneToCheck >= 0) && vecCore::Get( Active, laneToCheck );
    
    bool anyProgressLastStep = ! vecCore::MaskEmpty(goodStep) ;    
    cout << " - iter = " << iter << " anyProgressLastStep = " << anyProgressLastStep << std::endl;       
    ReportManyRowsOfDoubles("Current X/P",  yStart, 6 );
    ReportRowOfDoubles("Charge",  charge, 6 );       
    ReportManyRowsOfDoubles("d/ds [X,P] ", dydx, 6 );    
#endif    
    itersLeft--;
    iter++;
    if (partDebug) cout << " OneGoodStep - iteration = " << iter << "             "
                        << " ( total iterations = " <<  (tot_no_trials + iter) << " ) "
                        << endl;
    // #ifdef STORE_ONCE
    vecCore::MaskedAssign(h, finished, Real_v(0.0)); // Set h = 0.0 for finished lanes -- ensure no change !
                                                     // #endif

    GetStepper()->StepWithErrorEstimate(yStart, dydx, charge, h, ytemp, yerr); // CAREFUL -> changes for others ?
    //**********==>>>>>>>>>>>>>>>>>>>>>
    
    IncrementStepperCalls();     //  fStepperCalls++;

#ifdef DRIVER_PRINT_PROGRESS
    bool DebugEachIteration = false;
    if (partDebug && DebugEachIteration) {
      cout << "1st Report - after call to Step-With-Error-Estimate" << endl;
      // ReportRowOfBools<Real_v>( "(already) finished", finished );
      ReportRowOfBools<Real_v>("active", Active);
      // ReportManyRowsOfDoubles( "yStart",  yStart, Nvar );
      ReportRowOfDoubles("h", h);
      ReportManyRowsOfDoubles("yOut", ytemp, Nvar);
      ReportManyRowsOfDoubles("yerr", yerr, Nvar);
      // } else {
      // ReportRowOfDoubles( "h",          h );
    }
#endif

    ErrorEstimatorSixVec fErrorEstimator( eps_rel_max, GetMinimumStep() );
    Real_v magmom_sq = yStart[3] * yStart[3] + yStart[4] * yStart[4] + yStart[5] * yStart[5];
#ifdef CHECK_ONE_LANE
    Real_v errpos_sq = 0.0; // square of displacement error
    Real_v errmom_sq = 0.0; // square of momentum vector difference
    Real_v epsPosition = 0.0; 
    // errmax_sq = fErrorEstimator.EstimateError( yerr, h, magmomInit_sq, epsPosition,  errpos_sq, errmom_sq );
    errmax_sq = fErrorEstimator.EstimateError( yerr, h, magmom_sq, epsPosition,  errpos_sq, errmom_sq );
    if( 0 ) { // printLane ) {    
       std::cout << "Check after EstimateError:  errmom_sq= " << std::sqrt( vecCore::Get(errmom_sq, laneToCheck ));
       std::cout << " Current mom^2= " << vecCore::Get( magmom_sq, laneToCheck ) ;
       // std::cout << " Initial mom^2= " << vecCore::Get( magmomInit_sq, laneToCheck ) ;
       std::cout << "  yErr: ";
       for( unsigned int ii=0; ii<Nvar; ii++) 
          std::cout << "[" << ii << "] = " << vecCore::Get( yerr[ii], laneToCheck ) << " ";
       std::cout << std::endl;
    }
#else
    // errmax_sq = fErrorEstimator.EstimateError( yerr, h, magmomInit_sq );
    errmax_sq = fErrorEstimator.EstimateError( yerr, h, magmom_sq);
    // errmax_sq = fErrorEstimator.EstimateError( yerr, h, yStepEnd );    
#endif

#ifdef DRIVER_PRINT_PROGRESS    
    bool ReportIntegrationStep = false;
    if (partDebug && ReportIntegrationStep) {
      ReportRowOfDoubles("epsPositin", epsPosition);
      // ReportRowOfDoubles( "invEpsPos2", invEpsPositionSq );
      // ReportRowOfDoubles( "errpos_sq", errpos_sq );
      Real_v errPos = vecCore::math::Sqrt(errpos_sq);
      ReportRowOfDoubles("errPos", errPos);
      // ReportRowOfDoubles( "errmom_sq", errmom_sq );
      Real_v errMom = vecCore::math::Sqrt(errmom_sq);
      ReportRowOfDoubles("errMom", errMom);
      ReportRowOfSquareRoots("errmax", errmax_sq);
      ReportRowOfDoubles("errmax_sq", errmax_sq); // To compare with stored values
    }
#endif

    goodStep = Active && (errmax_sq <= 1.0);

#ifdef DRIVER_PRINT_PROGRESS        
    Real_v hAdvance= vecCore::Blend( goodStep, h, Real_v(0.0) );

    // Real_v xNext = x + hStep; //   Updates both good and bad steps -- else it should be  = x + hAdvance;
    Real_v xNext = x + hAdvance; //   Updates only good steps --  2019.01.07

    cout << " ====================================================================" << endl;
    ReportRowOfDoubles("charge ***>", charge );    
    ReportRowOfDoubles("x:   Start  =", x );
    ReportRowOfDoubles("xNext: End  =", xNext );
    ReportRowOfDoubles("hStep/tried =", h );    
    ReportRowOfDoubles("hAdvance/got=", hAdvance );
    // ReportRowOfDoubles("hdid        =", hdid );
    ReportRowOfBools<Real_v>("goodStep", goodStep);
    
    ReportRowOfDoubles("x: Updated =", x );    
    cout << " ====================================================================" << endl;
#endif 
    
    Bool_v laneDone = (goodStep | finished);
    bool allDone    = vecCore::MaskFull(laneDone);

    finished = laneDone;
    
#ifdef CHECK_ONE_LANE
    // Debugging one lane (at a time)  -------------  2019.02.27
    if( printLane ) { // if ( (laneToCheck >= 0) && vecCore::Get( Active, laneToCheck ) ) {
       ReportOneLane ( h, x,
                       epsPosition, errpos_sq, errmom_sq, errmax_sq, laneDone,
                       allDone, iter, tot_no_trials, laneToCheck, trackToPrint,
                       "SimpleID" );  // "SimpleIntDrv" );
       std::cout << " Track to check " << trackToPrint << " laneToCheck = " << laneToCheck << std::endl;       
       ReportRowOfSquareRoots("ErrPos", errpos_sq );
       ReportRowOfSquareRoots("ErrMom", errmom_sq );
       ReportRowOfSquareRoots("ErrMax", errmax_sq );
       ReportManyRowsOfDoubles("yErr", yerr, Nvar);
       
       if( 1 ) {
          std::cout << "SID: Status after stepper call ----------------------------------------------" << std::endl;
          FormattedReporter::FullReport(yStart, charge, dydx, h /*hStep*/, ytemp /*yStepEnd*/,
                                        yerr, errmax_sq, Active, goodStep );
       }
       if( 1 ) {
          // ReportRowOfSquareRoots("|err-p|", yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5] );
          // ReportRowOfDoubles("up = SumErr^2", sumerr_sq );
          // ReportRowOfDoubles("dwn= magMom^2+e", magmom_sq + tinyValue );
          // ReportRowOfDoubles("mul:1/e_vel^2", invEpsilonRelSq );
          // ReportRowOfDoubles("ErrMom^2", errmom_sq );
          // ReportRowOfSquareRoots("ErrMom", errmom_sq );
       }
    }
    // End debug code                 -------------  2019.02.27
#endif          

    Active   = !finished;
    
#ifdef DRIVER_PRINT_PROGRESS
    if (partDebug) {
      ReportRowOfBools<Real_v>("goodStep", goodStep);
      ReportRowOfBools<Real_v>("laneDone", laneDone);
      ReportRowOfBools<Real_v>("finished(upd)", finished);
      ReportRowOfBools<Real_v>("Active(upd)", Active);

      ReportRowOfDoubles("x: Updated =", x );
    }
#endif
    
    if (allDone) // All (or remaining) steps succeeded.
    {
      if (partDebug) cout << "SID: All done. Storing Good Values, then breaking." << endl;       
      // Idea 1.5
      StoreGoodValues(ytemp, h, errmax_sq, goodStep, yFinal, hFinal, errmax_sqFinal);
      //*************
      //  errmaxSqFallThru = errmax_sq;  
      break;
    }

    Real_v errPower = PowerSameIf<Real_v>(errmax_sq, fHalfPowerShrink, !laneDone);    
    Real_v hReduced = h * Max(Real_v(0.1), kSafetyFactor * errPower);

    Real_v hnew = vecCore::Blend(finished, Real_v(0.0), hReduced);
    xnew        = x + hnew;

    stepSizeUnderflow = Active && (xnew == x);

#ifdef CHECK_ONE_LANE    
    if( printLane ) {
      std::cout << "SID: After chance to break. (Not allDone.) Remaining lanes represent work: " << std::endl;
      ReportRowOfBools<Real_v>("laneDone", laneDone);
      ReportRowOfBools<Real_v>("Active", Active);
      ReportRowOfDoubles("powerShrink=", 2*fHalfPowerShrink );
      std::cout << " safetyFactor = " << fSafetyFactor << std::endl;
      ReportRowOfDoubles("errMaxSq=", errmax_sq );
      ReportRowOfSquareRoots("errMax=", errmax_sq );
      // ReportRowOfDoubles("errPower=", errPower );
      ReportRowOfDoublesIf<Real_v>("hReduced/ifActive", hReduced, Active);
      // ReportRowOfDoubles("hReduced=", hReduced );
      ReportRowOfDoubles("hnew=", hnew );
      ReportRowOfDoubles("xnew=", xnew );
      if( ! vecCore::MaskEmpty(stepSizeUnderflow) ) {
         ReportRowOfBools<Real_v>("underflow", stepSizeUnderflow);
         std::cout << "** WARNING: Underflow Detected in SimpleIntegrationDrive !!" << std::endl;
         std::cerr << "SimpleIntegrationDriver> Underflow Detected !!" << std::endl;
      }
    }
#endif

#ifndef STORE_ONCE
    if (!vecCore::MaskEmpty(stepSizeUnderflow)) {
      int numUnder = 0;
      for (unsigned int i = 0; i < vecCore::VectorSize<Real_v>(); i++) {
        if (vecCore::Get(stepSizeUnderflow, i)) {
          numUnder++;
        }
      }
      cout << "WARNING> Underflow detected in " << numUnder << " lanes." << endl;
      // Idea 1.0
      // StoreGoodValues( ytemp,      h,        errmax_sq,
      //                   stepSizeUnderflow, // && !finished,
      //                   yFinal,     hFinal,   errmax_sqFinal );
    }
    // Idea 1.5 :  Use only one store for goodStep & underflow lanes (if continuing)
    if (!vecCore::MaskEmpty(stepSizeUnderflow || goodStep)) {
#ifdef DRIVER_PRINT_PROGRESS       
      if (partDebug) cout << "Store and Report Stored lanes - v1.5 allDone - good or Underflow." << endl;
#endif
      //*************
      StoreGoodValues(ytemp, h, errmax_sq,
                      (goodStep || stepSizeUnderflow), // && !alreadyFinished,
                      yFinal, hFinal, errmax_sqFinal);
      //*************      
    }
#endif
    finished = finished || stepSizeUnderflow;

    h = hnew;

    // Interim Report of problem lanes -- For debugging only !!
    Bool_v problemLanes = stepSizeUnderflow; // && Active ;  // -> Already checking 'Active' above
    if (!vecCore::MaskEmpty(problemLanes)) {
      std::cerr << "GVIntegratorDriver::OneStep:" << std::endl
                << "  Stepsize underflow in Stepper ( Report 1 - in loop )" << std::endl;
      Bool_v problemLanes = stepSizeUnderflow && !finished;

      ReportConditionLanes(problemLanes, x, xnew, h, htry);
    }

#ifdef STORE_ONCE
    errmaxSqFallThru= errmax_sq;
#endif
    
  } while (itersLeft > 0 && (!vecCore::MaskFull(finished)) //  was MaskFull( stepSizeUnderflow || goodStep ) )
  );

  // CPU time for this next line grows much more than linearly with the number
  // of threads (eg from 0.0118% of total with 1 threads to 1.15% with 4 threads)
  // tot_no_trials += iter;

#ifdef STORE_ONCE
  //  'Idea 3' - Store exactly one time ( here - except on loop exit)
  StoreGoodValues(ytemp, h, errmaxSqFallThru, finished, yFinal, hFinal, errmax_sqFinal);
  //*************        
//   Why not store all ?
#endif

  if (!vecCore::MaskEmpty(stepSizeUnderflow)) {
    // int numUnder= NumberTrue( stepSizeUnderlow );
    std::cerr << "== Check after iteration loop: found " // << numUnder
              << "underflow lanes." << std::endl;
    ReportConditionLanes(stepSizeUnderflow, x, xnew, h, htry);
  }

#ifdef DRIVER_PRINT_PROGRESS  
  if (partDebug)
    cout << "SimpleIntDrv: 1-step - Loop done at iter = " << iter << " with h = " << h << " from htry= " << htry
         << std::endl;
#endif
  
  h = hFinal;

  // Check against fErrcon to avoid calling power ... saves work if any are 'over' max
  // fErrcon is 'const', cache them here to avoid memory contention (which leads to a lack of scalability)
  static const Real_v kErrCon2_v     = fErrcon * fErrcon;  
  static constexpr auto tPowerGrow_s = .5 * kPowerGrow;
  Bool_v underThresh                 = errmax_sqFinal <= kErrCon2_v;  
  Real_v errStretch1raw              = fSafetyFactor * PowerSameIf(errmax_sqFinal, tPowerGrow_s, !underThresh);
  // Real_v errStretch1raw           = fSafetyFactor * PowerSameIf(errmax_sqFinal, 0.5 * GetPowerGrow(), !underThresh);  

  // Note:  lanes with 'false' argument (i.e. underThresh=true) will have value 1.0
  Real_v errStretch =
    vecCore::Blend( underThresh, Real_v(fMaxSteppingIncrease), errStretch1raw );

  ReportRowOfDoubles("errMaxSq(Final)=", errmax_sqFinal );
  // ReportRowOfBools<Real_v>("underThresh", underThresh);
  std::cout << " errcon = " << errcon << std::endl;
  std::cout << " Power-Grow = " << 0.5 * GetPowerGrow() << std::endl;
  if( errcon == 0.0 ) { std::cerr << " Simple IntDrv> ERROR: errcon is Zero.  Value = " << errcon << std::endl; } 
  
  hnext = errStretch * h;

#define CHECK_STRETCH_FACTOR 1  
#ifdef CHECK_STRETCH_FACTOR
  constexpr double minErr2 = 1e-100;
  Real_v emax2pos          = Max(errmax_sqFinal, Real_v(minErr2));
  // Real_v errStretchOpen     = fSafetyFactor * Exp((0.5 * GetPowerGrow()) * Log(emax2pos)); // Was: Log(errmax_sqFinal) );
  Real_v errStretchOpen     = fSafetyFactor * Math::Pow( emax2pos , Real_v(0.5 * GetPowerGrow()) );
  // ReportRowOfDoubles("stretch-raw/new", errStretch1raw);  
  // ReportRowOfDoubles("stretch-raw/old", errStretchOpen);
  Real_v errStretchOld  = Min( errStretchOpen , Real_v(fMaxSteppingIncrease) );
  // ReportRowOfDoubles( "old: errStretch (initial) ", errStretchOld);
  Bool_v zeroErr = errmax_sqFinal <= minErr2;
  vecCore::MaskedAssign(errStretchOld, zeroErr, Real_v(fMaxSteppingIncrease));
  // ReportRowOfDoubles( "old: errStr (change4tiny) ", errStretchOld);

  if (!vecCore::MaskEmpty( vecCore::math::Abs(errStretch - errStretchOld) > 1.0e-9 * errStretch))
  {
    cout << "################################################################################-----------------" << std::endl;
    cout << "ERROR> Lanes found with differences in calculated value of 'errStretch'"
         << "       ( factor for stretching step size for 'good' next step." << endl;
    ReportRowOfDoubles("new-old: errStretch", errStretch - errStretchOld);
    ReportRowOfDoubles("old: errStretch", errStretchOld);
    ReportRowOfDoubles("new: errStretch", errStretch);
    // ReportRowOfDoubles("errStretch-raw", errStretch1raw);
    cout << "################################################################################-----------------" << std::endl;    
  }
#endif

#ifdef CHECK_ONE_LANE
  // bool  printDbgH = (laneToCheck >= 0) && vecCore::Get( htry, laneToCheck ) > 0.0 ;
  bool  printDbgH = false;  // (laneToCheck >= 0) && vecCore::Get( htry, laneToCheck ) > 0.0 ;    
  if( printDbgH ) {
      std::cout << "################################################################################-----------------" << std::endl;
      std::cout << "Determining next step size (hnext) as steps were successful (in all lanes.)" << std::endl;
      ReportRowOfDoubles("errMaxSq(Final)=", errmax_sqFinal );
      // ReportRowOfDoubles("errStretch-raw", errStretch1raw);
      ReportRowOfDoubles("errStretch", errStretch);
      // ReportRowOfDoubles("old: errStretch", errStretchOld);
      ReportRowOfDoubles("hnext=", hnext );      
      std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;
      ReportRowOfSquareRoots("sqrt(errMax2/fErrcon^2)=", (1.0 / (errcon * errcon)) *  errmax_sqFinal);
      ReportRowOfDoubles(    "errMax/Errcon (by hand)=", vecCore::math::Sqrt( errmax_sqFinal ) / errcon );
      // ReportRowOfBools<Real_v>("overThresh", overThresh);
      ReportRowOfBools<Real_v>("underThresh", underThresh);
      std::cout << "***********   errcon = " << errcon << "    square = " << errcon * errcon << std::endl;
      ReportRowOfDoubles("powerGrow=", GetPowerGrow() );
      std::cout << " safetyFactor = " << fSafetyFactor << std::endl;
      std::cout << "################################################################################------------------" << std::endl;      
  }
#endif
  
  // ReportRowOfDoubles( "OGS: h-final", hFinal);

  hdid = hFinal;
  x += hdid;

#ifdef DRIVER_PRINT_PROGRESS  
  bool OGSreport = true;
  if (partDebug && OGSreport) {
    ReportRowOfDoubles("x: Updated =", x );
     
    ReportRowOfDoubles("OGS: errmax2", errmax_sqFinal);
    ReportRowOfSquareRoots("OGS: errmax    ", errmax_sqFinal);
    ReportRowOfDoubles("OGS: x (new)   ", x);
    ReportRowOfDoubles("OGS: h-did     ", hdid);
    ReportRowOfDoubles("OGS: h-next    ", hnext);
    ReportRowOfDoubles("OGS: facStretch", errStretch);

    Real_v xFinish = xStart + hdid;    
    ReportRowOfDoubles("KPS: x0+hdid ", xFinish);
    // ReportRowOfDoubles( "OGS: hFinal", hFinal);
  }

  if (partDebug) {
    cout << " hdid= " << hdid << " and hnext= " << hnext << std::endl
         << " end of  OneGoodStep method .......................... " << std::endl;
  }
#endif
  // for(int k=0;k<Nvar ;k++) { y[k] = yFinal[k]; }
  
  return;
} // end of OneGoodStep()

// ---------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
inline void SimpleIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::StoreGoodValues(
    const Real_v yWork[], const Real_v &hValue, const Real_v &errMaxSqValue, const vecCore::Mask_v<Real_v> &storeFlag,
    Real_v yFinal[], Real_v &hFinalStore, Real_v &errMaxSqStore) const
// yWork,  hValue,      epsSqVal   represent the output variables
// yFinal, hFinalStore, epsSqStore represent the output variables
{
  using std::cout;
  using std::endl;
  using vecCore::Get;
  using vecCore::MaskedAssign;
  using vecCore::Set;

  if (vecCore::MaskFull(storeFlag)) {
    // std::cout << "StoreGoodValues: Unconditional assignment to output - all together." << std::endl;

    for (unsigned int j = 0; j < Nvar; ++j)
      yFinal[j] = yWork[j];

    hFinalStore   = hValue;
    errMaxSqStore = errMaxSqValue;
  } else {
    for (unsigned int j = 0; j < Nvar; ++j)
      MaskedAssign(yFinal[j], storeFlag, yWork[j]);

    MaskedAssign(hFinalStore, storeFlag, hValue);
    MaskedAssign(errMaxSqStore, storeFlag, errMaxSqValue);
  }
  // All the work has been done now.

#ifdef DRIVER_PRINT_PROGRESS  
  // Print the input & output
  bool verboseStore = false;
  if (verboseStore && partDebug) {
    using FormattedReporter::ReportManyRowsOfDoubles;
    using FormattedReporter::ReportRowOfBools;
    using FormattedReporter::ReportRowOfDoubles;
    cout << "==============================================" << endl;
    cout << "Called Store-Final-Values.  Input is " << endl;
    ReportRowOfDoubles("h", hValue);
    ReportManyRowsOfDoubles("y", yWork, Nvar);
    ReportRowOfDoubles("errmaxSq", errMaxSqValue);

    ReportRowOfBools<Real_v>("storeFlag", storeFlag);
    cout << "==============================================" << endl;
    cout << "Report of Stored lanes - in Store-Final-Values" << endl;
    // ReportArray( methodName, "hFinal", hFinal );
    ReportRowOfDoubles("hFinal", hFinalStore);
    ReportManyRowsOfDoubles("yFinal", yFinal, Nvar);
    ReportRowOfDoubles("errmaxSq/Final", errMaxSqStore);
    // cout<< " yerr is: " << yerr[0] << endl;
  }
#endif
  
#ifdef CHECK_STORE
  bool checkFlag = true;
  if (checkFlag) {
    for (unsigned int i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
      if (storeFlag[i]) {
        if (partDebug)
          cout << "Checking store of lane " << i << " :  "
               << " h = " << Get(hValue, i)               // hValue[i]
               << " errMaxSq = " << Get(errMaxSqValue, i) // errMaxSqValue[i]
               << " ... ";                                // << endl;
        //  Check  hFinalStore [i] =      hValue[i];
        double hStored      = Get(hFinalStore, i);
        double hOriginal    = Get(hValue, i);
        const double epsTol = 1.0e-9;
        assert(std::fabs(hStored - hOriginal) <= epsTol * std::fabs(hOriginal));
        //  Check errMaxSqStore [i] == errMaxSqValue[i];
        double emsStored   = Get(errMaxSqStore, i);
        double emsOriginal = Get(errMaxSqValue, i);
        assert(std::fabs(emsStored - emsOriginal) <= epsTol * std::fabs(emsOriginal));
        for (int j = 0; j < ncompSVEC; ++j) {
          double yStored   = Get(yFinal[j], i);
          double yOriginal = Get(yWork[j], i);
          assert(std::fabs(yStored - yOriginal) <= epsTol * std::fabs(yOriginal));
        }
        if (partDebug) cout << " .. OK " << endl;
      } // else // Check that they are unchanged ... ?
    }
  }
#endif

#ifdef ORIGINAL_CODE
  for (int i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
    if (storeFlag[i]) {
      cout << "Storing lane " << i << " :  "
           << " h = " << Get(hValue, i)                        // hValue[i]
           << " errMaxSq = " << Get(errMaxSqValue, i) << endl; // errMaxSqValue[i]

      // finishedArr[i] = true;
      // Set( finished, i, true);  //  finished[i] = true;

      //   hFinalStore [i] =      hValue[i];
      Set(hFinalStore, i, Get(hValue, i));
      //   errMaxSqStore [i] =    errMaxSqValue[i];
      Set(errMaxSqStore, i, Get(errMaxSqValue, i));
      for (int j = 0; j < ncompSVEC; ++j) {
        //   yFinal[j] [i] =     yWork[j] [i];
        Set(yFinal[j], i, Get(yWork[j], i));
      }
    }
  }
#endif
}


// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
int SimpleIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::InitializeLanes(
    const FieldTrack yInput[],
    const double     hstep[],
    const double     charge[],
    int              nTracks,
    int       indexArr[], // [vecCore::VectorSize<Real_v>()]
    Real_v    y[],
    Real_v  & hStepLane,
    Real_v  & chargeLane,
    Real_v  & startSlen,
    int     & numFilled, // How many were loaded.
    bool      badStepSize[],  // Output
    int     & numBadSize ) const
// Converts input scalar stream to acceptable form of Vc vectors
// for vector processing in OneStep
{
  using FormattedReporter::ReportArray;
  using vecCore::Set;
  // void ReportArray( const std::string& context, const std::string& varName,
  //                  const double Arr[],         int numTracks,             bool banner= false );

  if (partDebug) std::cout << "----Initializing Lanes ----" << std::endl;
  // const int NumComptFT= FieldTrack::NumCompFT;
  double yStartScalar[Nvar]; // [ncompSVEC]; // Size depends on needs of DumpToArray

  if (partDebug) {
    std::cout << "InitLanes:  chargeArr [0] = " << charge[0] << " ,  [1] = " << charge[1] << ",  [2] = " << charge[2]
              << ",  [3] = " << charge[3] << std::endl;
    ReportArray("InitLanes", "hStepArr", hstep, nTracks, true);
    ReportArray("InitLanes", "chargeArr", charge, nTracks);
  }
  hStepLane  = Real_v(-12345.6789); //  Signals 'inactive' for lanes not loaded
  chargeLane = Real_v(0.0);         //    >>  Ditto  >>

  // Rest only at the start 
  numBadSize  = 0; // Ensure it starts at Zero !?
  int j       = 0;
  
  size_t slot = 0;
  do {
    double hVal       = hstep[j];
    double chargeVal  = charge[j];
    bool invalidTrack = (hVal <= 0.0) || (chargeVal == 0.0);
    badStepSize[j] =  invalidTrack;

    if (invalidTrack) {
      if (hVal <= 0.0)
        std::cout << " WARNING> Non-positive step-size h = " << std::setw(5) << hVal << " for j= " << j << std::endl;
      if (chargeVal == 0.0) std::cout << " WARNING> Zero charge " << chargeVal << " for j= " << j << std::endl;
      // Will need to copy input to output ...
      numBadSize++;
    } else {
      indexArr[j] = slot;
      yInput[j].DumpToArray(yStartScalar);

      //   hStepLane   []  = hVal;
      Set(hStepLane, slot, hVal);
      //   chargeLane  []  = chargeVal ;  // -- i.e. charge[j];
      Set(chargeLane, slot, chargeVal); // charge[j] );
      // std::cout << " slot = " << slot << " j= " << j << " charge = " << charge[j] << std::endl;
      //  startSlen  [j] = xStart[j] ); // yInput[j].GetCurveLength();
      Set(startSlen, slot, // xStart[j] );
          yInput[j].GetCurveLength());

      for (unsigned int i = 0; i < Nvar /*ncompSVEC*/ ; ++i) {
        //   y[i] [j] = yStartScalar[i];
        Set(y[i], slot, yStartScalar[i]);
      }

#ifdef CHECK_ONE_LANE
      const int trackToPrint = IntegrationDriverConstants::GetInstance()->GetTrackToCheck();      
      // if( j == trackToPrint ) { laneToCheck = slot; }
      if( j == trackToPrint ) {
         // if( laneToCheck == -1 ) ...
         std::cout << "SID::InitLanes> found lane = " << slot
                   << " (was = " << laneToCheck << " ) "
                   << " for trackToPrint = " << trackToPrint << std::endl;
         laneToCheck = slot;
         // std::cout << "SID::InitLanes> found lane = " << laneToCheck << " for trackToPrint = " << trackToPrint << std::endl;         
      }
#endif
      ++slot;
    }
    ++j;
  } while (slot < vecCore::VectorSize<Real_v>() && j < nTracks);

  numFilled = slot;
  return j; // Next location, after all used to fill

} // End of InitializeLanes function

// ------------------------------------------------------------

template <typename Real_v, typename scalar_t>
bool CheckLaneValue(Real_v varLanes, int lane, scalar_t expected, std::string testName, std::string nameVal)
{
  double current = vecCore::Get(varLanes, lane);
  bool ok        = (current == expected);
  if (!ok) {
    std::cerr << testName << " : ERROR in comparison of " << nameVal << " [ lane = " << lane << " ] "
              << " current = " << current << " VS expected = " << expected << "  diff = " << current - expected
              << std::endl;
  }
  return ok;
}

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
bool SimpleIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::InsertNewTrack(
    const FieldTrack yInput[], const double hstep[], const double charge[], const int slot, int &trackNextInput,
    bool       succeeded[],      // Output 
    Real_v     y[], // [Nvar]
    Real_v   & hStepLane,
    Real_v   & chargeLane,
    Real_v   & startCurveLength,
    int        indexArr[], // [vecCore::VectorSize<Real_v>()]
    int        nTracks,
    bool       badStepSize[],  
    int      & numBadSize
   ) const

// Inserts a new track whenever a lane is finished.
// returns isDoneLane = true for h<=0 case, false otherwise
// because in former case, no further work is required
{
  using vecCore::Set;

  if (partDebug) std::cout << "----Inserting New Track " << trackNextInput << " at position " << slot << std::endl;

  bool filled = false; // to get the while loop starting
  while (trackNextInput < nTracks && !filled) {
    // Ensure that hstep > 0
    double hStepNext = hstep[trackNextInput];

    bool goodStep = (hStepNext > 0.0) && (charge[trackNextInput] != 0.0);
    badStepSize[trackNextInput] =  !goodStep;
    
    if ( goodStep ) {
      double yScalar[Nvar /* fNoVars*/ ];

      yInput[trackNextInput].DumpToArray(yScalar);
      // for (int i = 0; i < Nvar; ++i)
      // const int NumComptFT= FieldTrack::NumCompFT;
      for (unsigned int i = 0; i < Nvar /*ncompSVEC*/ ; ++i) {
        //   y[i] [slot] = yScalar[i];
        Set(y[i], slot, yScalar[i]);
      }
      indexArr[slot] = trackNextInput;
      //   hStepLane  [slot] = hstep[trackNextInput];
      Set(hStepLane, slot, hstep[trackNextInput]);
      //   chargeLane [slot] = charge[trackNextInput] );
      Set(chargeLane, slot, charge[trackNextInput]);

      double slen = yInput[trackNextInput].GetCurveLength();
      //    startCurveLength [slot] = slen;
      Set(startCurveLength, slot, slen);

#ifdef CHECK_ONE_LANE
      const int trackToPrint = IntegrationDriverConstants::GetInstance()->GetTrackToCheck();
      if( trackNextInput == trackToPrint )
         laneToCheck = slot;
#endif

      filled = true;
    } else {
      // A zero or negative step is anomalous - an error
      succeeded[trackNextInput] = (hStepNext == 0.0);
      numBadSize++;
      
      if (hStepNext == 0) {
         std::cerr << "Proposed step is zero; hstep = " << hStepNext;
      } else {
        std::cerr << "Invalid run condition." << std::endl
                  << "Proposed step is negative; hstep = " << hStepNext << ".";
      }
      std::cerr << " for track " << trackNextInput
                << " with position/momentum " << yInput[trackNextInput]
                << std::endl;
    }
    trackNextInput++;
  }

  return filled;

} // End of InsertNewTrack function

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void SimpleIntegrationDriver< T_Stepper, Nvar>::
   StoreOutput(const Real_v  yEnd[],
               const Real_v  x,
               int           currIndex,
               FieldTrack    yOutput[],
               int           indOut,
               const double  hstep[],
               bool          succeeded[],
               int           nTracks
   ) const
// Called whenever a lane is finished, to store output of integration
// Input arguments
//    - currIndex is the index of finished lane in Vc vector
//    - yEnd[] vector of values of integrands
//    - x      vector of final length (independent parameter value)
//    - hstep  don't store if h<=0
// Stores
//    - end position and momentum in yOutput ('scalar' output array)
//    - final curve length,        and
//    - success flag in array succeeded[nTracks]
//
{
  if (partDebug)
    std::cout << "----Storage position (out-arr): "
              << indOut
              // << " (ntracks= " << nTracks << ")"
              << std::endl;

  (void)nTracks; // Use the value in case on non-Debug builds - avoids compilation warning

  assert(0 <= indOut && indOut < nTracks && "Track Index is Out of Range");
  assert(0 <= currIndex && ((unsigned long)currIndex < vecCore::VectorSize<Real_v>()) && "Lane Index is Out of Range");

  double hOriginal = hstep[indOut];

  if (hOriginal >= 0.0) {
    // need to get a yEnd : scalar array
    double yOutOneArr[Nvar]; // std::max(int(ncompSVEC ), int(Nvar))];
    for (unsigned int i = 0; i < Nvar; ++i) // Was: i < fNoIntegrationVariables; ++i)
    {
      // yOutOneArr[i] =    yEnd[i] [currIndex]; // Constant col no., varying row no. for required traversal
      yOutOneArr[i] = vecCore::Get( yEnd[i], currIndex);
    }
    yOutput[indOut].LoadFromArray(yOutOneArr, Nvar);            // fNoIntegrationVariables);
    yOutput[indOut].SetCurveLength(vecCore::Get(x, currIndex)); // x is a double_v variable
  } else {
    succeeded[indOut] = (hOriginal == 0.0);
  }

#ifdef CHECK_ONE_LANE
  const int trackToPrint = IntegrationDriverConstants::GetInstance()->GetTrackToCheck();              
  if( indOut == trackToPrint )
     laneToCheck = -1;
#endif
  
} // End of StoreOutput function

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void SimpleIntegrationDriver<T_Stepper, Nvar>::AccurateAdvance(const FieldTrack yInput[],
                                                               const double     hstep[],
                                                               const double     charge[],
                                                               // double           epsilon, // Can be scalar or varying
                                                               FieldTrack       yOutput[],
                                                               bool             stillOK[],
                                                               int              nTracks) const
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

  static constexpr const char *methodName = "SID::AccurateAdvance";
  using FormattedReporter::GetMomentumMag;
  using FormattedReporter::ReportArray;
  using FormattedReporter::ReportManyRowsOfDoubles;
  using FormattedReporter::ReportRowOfBools;
  using FormattedReporter::ReportRowOfDoubles;
  using FormattedReporter::ReportRowOfSquareRoots;
  using FormattedReporter::ReportRowsOfPositionsMomenta;
  using FormattedReporter::GetMomentumMag;
  using FormattedReporter::ReportArray;
// #endif
  using Bool_v = vecCore::Mask_v<Real_v>;

  // using AllignedInt_v = vecCore::Int_v<Real_v>; // Alligned with Real_v ... same # of entries ..
  using vecCore::math::Min;
  using vecCore::math::Max;
  using std::cout;
  using std::endl;

  const double epsilon= FlexIntegrationDriver::GetMaxRelativeEpsilon();
  
  constexpr unsigned int VecSize = vecCore::VectorSize<Real_v>();
  int indexArr[VecSize]; // vecCore::VectorSize<Real_v>()];

  // Working variables for integration
  Real_v x, hnext, hdid, chargeLane, xStartLane(0.0), hStepLane;
  Real_v y[Nvar];     // Working array 1
  Real_v yNext[Nvar], // Working array 2
         dydx[Nvar];

  // Real_v ySubStepStart[Nvar /*ncompSVEC*/ ];           
  bool badStepSize[nTracks];    // Step size is trivial (0) or 'illegal' ( h < 0 )

  std::fill_n( stillOK,      nTracks, true); 
  // std::fill_n(stillOK, nTracks, 1); //  Assume success, flag failure ... Strange/JA
  // std::fill_n( badStepSize,  nTracks, false);
#ifdef DRIVER_PRINT_PROGRESS    
  if (partDebug) ReportArray(methodName, "hstep", hstep, nTracks);
#endif

  Bool_v /*lastStepOK,*/ succeededLane(false), isLastStepLane(false);
  Bool_v isDoneLane(false); // set true when there is a return statement
  int numFilled  = 0;
  int numBadSize = 0; // How many tracks have invalid step size
  // ThreadLocal int  noGoodSteps =0 ;  // Bad = chord > curve-len

  for (unsigned int i = 0; i < VecSize; ++i)
    indexArr[i] = -1;

  int idNext = InitializeLanes(yInput, hstep, charge, nTracks,                   // Input
                               indexArr,                                         // out - indices of lanes
                               y, hStepLane, chargeLane, xStartLane, numFilled,  // Main output - values
                               badStepSize, numBadSize);                         // 'Bad' output
  // OLD:
  // int     idNext = vecCore::VectorSize<Real_v>(); // ASSUMES that nTracks >= vecCore::VectorSize<Real_v>() !!! FIX

  if( numBadSize > 0 ) 
     cout << " InitializeLanes returned with " << numBadSize << " 'bad' tracks - i.e. hStep <= 0" << endl;  

  if( numBadSize > 0)
   for (int i = 0; i < idNext; ++i)
      if( badStepSize[i] )
      {
         yOutput[i] = yInput[i];
      }

  
  // const maxSize = std::max( (int) vecCore::VectorSize<Real_v>(), (int) numTracks );
  if (idNext > (int)vecCore::VectorSize<Real_v>()) // Some lanes had hstep <= 0.0
  {
    for (unsigned int i = 0; i < VecSize; ++i) // vecCore::VectorSize<Real_v>(); ++i)
    {
      if (hstep[i] <= 0.0) {
        stillOK[i] = (hstep[i] == 0.0);
      }
    }
  }

  succeededLane = (hStepLane <= 0.0);

  // Assume that all lanes are currently Full ?? ===>>> TEST it !!

  using Index_v = vecCore::Index<Real_v>;
  Index_v nstp(0); // Should be Int_v with size compatible with Real_v ( likely 64-bit )

#ifdef DRIVER_PRINT_PROGRESS       
  // Real_v yStart[Nvar];    // Starting values
  // for( int i=0; i<Nvar; i++ ) { yStart[i] = y[i]; }
  // Real_v momStart = GetMomentumMag(yStart);  // Can be done later, and if needed
  Real_v momStart = GetMomentumMag(y);  // True starting momentum
#endif
  
  // Real_v StartPosAr[3];
  // StartPosAr[0] = y[0];   StartPosAr[1] = y[1];  StartPosAr[2] = y[2];

  // isDoneLane needed. In end, other conditions might keep changing
  // even if processing for that lane is finished.
  // Need a way to store the fact that the lane is done.
  // Either make a new single condition that combines isDoneLane
  // and all other conditions or some conditions at least
  // For now, just adding isDoneLane : needs to be && or || with the first 3
  // and keep nTracks condition in final ||
  // Say continue if isDoneLane is not 1111 and rest all conditions are not 0000
  // while ( !vecgeom::IsEmpty((nstp<=fMaxNoSteps) && (x < x2) && (!isLastStepLane)) || idNext < nTracks  )

  Real_v hTry  = hStepLane;
  Real_v x1 = xStartLane;
  Real_v x2 = x1 + hStepLane;   // End of integration interval - value renewed as needed in refill
  
  x = x1;

  // const AllignedInt_v laneNum =

  auto pStepper = GetStepper();
  // TStepper_t * Stepper = GetStepper();  

  while ( ( !vecCore::MaskFull(isDoneLane)
            && !vecCore::MaskEmpty((nstp <= GetMaxNoSteps() )
                                   && (x < x2) && (!isLastStepLane)))
          ||
             idNext < nTracks
     )
  {     
#ifdef DRIVER_PRINT_PROGRESS     
    if (partDebug) {       
      std::cout << "******************************************************" << std::endl
                << "** Top of while loop (AccurateAdvance) ***************" << endl
                << "----hStepLane is: " << hStepLane << endl;
      //  ReportManyRowsOfDoubles("yStart", y, Nvar);
      // ReportManyRowsOfInts("Lane#", laneNum, Nvar);      
    }
#endif

    // if( h > fMinimumStep ) { QuickAdvance .. } else { .. below  //  ( Sequential code  )
    // if (useOneStep) {
    // GetStepper()
    pStepper ->RightHandSideInl( y, chargeLane, dydx ); // TODO: change to inline
    //---------****************-----------------------

#ifdef DRIVER_PRINT_PROGRESS         
    // using FormattedReporter::ReportManyRowsOfDoubles;
    // using FormattedReporter::ReportRowOfDoubles;
    ReportRowOfDoubles("hTry(before)", hTry);
    
    if( 0 ) {
       cout << "##-------------------------------------------------------------------------------" << endl;           
       cout << "### Accurate Advance ---- After RightHandSide" << endl;
       /* FormatterReported::*/ ReportManyRowsOfDoubles /*<Real_v>*/ ("dydx", dydx, Nvar);
       /* FormatterReported::*/ ReportRowOfDoubles("charge", chargeLane);
       cout << "##-------------------------------------------------------------------------------" << endl;       
    }
#endif
    
    OneGoodStep<Real_v>(y, dydx, chargeLane, x, hTry, epsilon, yNext, hdid, hnext);
    //*********---------------------------------------------------------
    // } else KeepStepping( y, dydx, x, h, epsilon, hdid, hnext, hStepLane, hDone) ;

#ifdef DRIVER_PRINT_PROGRESS    
    if (partDebug) {
      cout << "##-------------------------------------------------------------------------------" << endl;
      cout << "### Accurate Advance ---- After return from OneGood Step" << endl;
      ReportRowOfDoubles("charge", chargeLane);
      ReportRowOfDoubles("hTry (after)", hTry);      
      // ReportRowOfDoubles("h-ask", h);
      ReportRowOfDoubles("h-did", hdid);
      ReportRowOfDoubles("hNext", hnext);
      cout << "##-------------------------------------------------------------------------------" << endl;
      ReportRowOfDoubles("x", x);
      ReportRowOfDoubles("xExpectEnd/x2", x2);
      cout << "##-------------------------------------------------------------------------------" << endl;
      // ReportManyRowsOfDoubles("yStart/True", yStart, Nvar);
      ReportManyRowsOfDoubles("yStepStart", y, Nvar);
      ReportManyRowsOfDoubles("dydx", dydx, Nvar);
      cout << "##-------------------------------------------------------------------------------" << endl;
      // ReportManyRowsOfDoubles("yNext", yNext, Nvar);
      // Real_v momStart = GetMomentumMag(yStart); // True initial momentum magnitude      
      Real_v momStepStart = GetMomentumMag(y);  //  'Step' start - easier to have
      ReportRowsOfPositionsMomenta("yNext", yNext, Nvar, momStepStart);
      // Real_v momTrueStart = GetMomentumMag(yStart); // True initial momentum magnitude
      // ReportRowsOfPositionsMomenta("yNext", yNext, Nvar, momTrueStart);  // Check vs true initial mag.
      // cout << "##-------------------------------------------------------------------------------" << endl;      
      // ReportManyRowsOfDoubles( "yNext", yNext, Nvar);
    }
#endif

    // lastStepOK = (hdid == h);
    // ++fNoTotalSteps;

#ifdef DRIVER_PRINT_PROGRESS
    static constexpr bool reportMove = true;    
    if (partDebug && reportMove) {
      ReportRowOfDoubles("charge", chargeLane);
      
      // ThreeVector EndPos( y[0], y[1], y[2] ); // Check the endpoint
      const Real_v edx = yNext[0] - y[0], edy = yNext[1] - y[1], edz = yNext[2] - y[2];
      Real_v endPointDist = vecgeom::Sqrt(edx * edx + edy * edy + edz * edz);
      ReportRowOfDoubles("Move-x", edx);
      ReportRowOfDoubles("Move-y", edy);
      ReportRowOfDoubles("Move-z", edz);
      ReportRowOfDoubles("Move-L", endPointDist);
      ReportRowOfDoubles("Move-L/hdid", endPointDist / hdid);
    }
#endif

    // Note: xStartLane must be positive. ( Ok - expect starting value = 0 )
    Real_v stepThreshold           = vecCore::math::Min(epsilon * hStepLane, GetSmallestFraction() * (xStartLane + hdid));
    Bool_v avoidNumerousSmallSteps = hTry < stepThreshold;

    // If it is always true for h<=0 --> lastStep is true, hence the lane will be sent to StoreOutput.

    isLastStepLane = avoidNumerousSmallSteps || isLastStepLane;
    // 'Accumulate' ie use 'OR' - in case lane is already finished or empty

    // x += hdid;  It is already updated - do not add the step again!!

    x2              = x1 + hStepLane;
    Real_v xremains = x2 - x; // (hStepLane - x) + x1; // was ( x2 - x )
    // For rest, check the proposed next stepsize

#ifdef DRIVER_PRINT_PROGRESS
    if (partDebug) {
      cout << " hRequest= " << hStepLane << endl;
      // cout << " x-Start = " << xStartLane << endl;
      cout << " hdid    = " << hdid << endl;
      cout << " x (Now) = " << x << endl;
      cout << " x2 -x   = " << xremains << endl;
    }
#endif

    hnext = Max(hnext, Real_v(GetMinimumStep()));
    // Note: This has potential for instability i.e. if MinStep is 'too long'

    // Ensure that the next step does not overshoot
    hnext = Min(xremains, hnext);

#ifdef DRIVER_PRINT_PROGRESS
    if ( 1 /*partDebug*/ ) {
       cout << "AccurateAdvance: hnext = " << hnext << " to replace hTry = " << hTry << endl;       
       ReportRowOfDoubles("RID/AccAdv: hNext=", hnext);
    }    
#endif

    hTry = hnext;

    // When stepsize overshoots, decrease it!
    // Must cope with difficult rounding-error issues if hstep << x2

    isLastStepLane = (hTry == 0.0) || isLastStepLane;

    if (partDebug) std::cout << " lastStep : " << isLastStepLane << std::endl;

    nstp += 1; // nstp++;

    succeededLane = (xremains <= 0.0); // (x>=x2); // If it was a "forced" last step ?

    Bool_v laneContinues = (nstp <= GetMaxNoSteps() ) && !succeededLane && !isLastStepLane;

    Bool_v renewedLanes(false); // To be 'set' only in the slots in which new values are inserted

    // Prepare values for next loop step.  Needed only if some lanes are continuing
    if (!vecCore::MaskEmpty(laneContinues)) // At least one lane continues
    {
      for (unsigned int i = 0; i < Nvar; ++i) {
        y[i] = yNext[i];
      }
      if (partDebug) cout << "Copying      y <- yNext , as 1+ lanes continue." << endl;
    } else {
      if (partDebug) cout << "Did NOT copy y <- yNext , as NO lanes continue." << endl;
    }

    if (!vecCore::MaskFull(laneContinues)) // At least one lane is finished
    {
      Bool_v finishedLane = !laneContinues;

#ifdef DRIVER_PRINT_PROGRESS
      if (partDebug) {
        cout << "SiD: At least one lane finished " << std::endl;
        cout << "  finishedLane        : " << finishedLane << std::endl;
        Bool_v CondNoOfSteps = (nstp <= Index_v(fMaxNoSteps));
        cout << "  Cond numSteps < Max : " << /* (nstp<=fMaxNoSteps)*/ CondNoOfSteps << std::endl;
        cout << "  Cond    (x < x2)    : " << !succeededLane << std::endl;
        cout << "  Cond  not Last Step : " << !isLastStepLane << std::endl;
      }
#endif

      // #if SINGLE_INSERT
      // Use a single method (to be written) to insert all new tracks ... avoid multiple 'calls' ?
      // InsertSeveralTracks( yInput, hstep,     charge,     i, idNext, stillOK,
      //                      y,      hStepLane, chargeLane, xStartLane );  .......
      // #else // SINGLE_INSERT
      for (unsigned int i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
        using vecCore::Set;
        if (partDebug)
          cout << "Checking for storing lane [ " << i << " ]  nstp = " << vecCore::Get(nstp, i)
               << " <= ? ( = fMaxNoSteps ) " << endl;

        // if ( finishedLane[i] &&  indexArr[i] != -1)
        if (vecCore::Get(finishedLane, i) && indexArr[i] != -1)
        {
          // 1. Store the Results (Output)
          stillOK[indexArr[i]] = vecCore::Get(succeededLane, i); // succeededLane[i];
          if (partDebug) std::cout << "----Storing Output at position: " << i << std::endl;
          StoreOutput(yNext, x, i, yOutput, indexArr[i], hstep, stillOK, nTracks); // Second - can change 'succeeded'
          //*********----------------------------------------------
          // Ananya: Do not pass succeededLane to StoreOutput (preference?), so
          //         'stillOK' should *not* be absorbed in StoreOutput.

          // 2. Load more work (if some exists) into empty slots
          //    TODO-1:  if all lanes are empty, can load in 'Vector mode' (aligned together)
          if (idNext < nTracks) {
            bool filled = InsertNewTrack(yInput, hstep, charge, i, idNext, stillOK, // Input values
                                         y,  hStepLane, chargeLane, xStartLane,    // Output values
                                         indexArr, nTracks,
                                         badStepSize,  numBadSize );

            Set(   isDoneLane, i, !filled); // isDoneLane[i]   = !filled;
            Set( finishedLane, i, !filled); // finishedLane[i] = !filled;
            Set( renewedLanes, i,  filled);
            // Set(     momStart, i, GetMomentum(y[i] )  );  // i= laneUsed

#ifdef DRIVER_PRINT_PROGRESS
            // for( unsigned int j=0; j<Nvar; j++ ) { Set(  yStart[j], i, y[j] ); }
            
            bool reportInsertNewTrack = true;
            if (partDebug && reportInsertNewTrack) {
              cout << " --Inserting New Track - part 1/2: loaded new state: " << (filled ? " Yes " : " No  ");
              cout << endl;
              ReportRowOfBools<Real_v>("renewedLanes", renewedLanes);
              ReportRowOfDoubles("hstep", hStepLane);
              ReportManyRowsOfDoubles("yCurrent", y, Nvar);
              ReportRowOfDoubles("xCurrent", x);
              // ReportRowOfDoubles( "xStart",  xStartLane);
              ReportRowOfDoubles("charge", chargeLane);
            }
#endif            
          } else {
            //   isDoneLane [i] = true;
            Set(isDoneLane, i, true);
            //   renewedLanes [i] = false;      // Not needed - it's the starting value
            // Set( renewedLanes, i,   false);
            indexArr[i] = -1;
            if (partDebug) cout << " --No New Tracks available: idNext = " << idNext << " nTracks= " << nTracks << endl;
          }
        }
      } // for ( uint i = 0; i < vecCore::VectorSize<Real_v>(); ++i)
      // #endif  // SINGLE_INSERT

      if (!vecCore::MaskEmpty(renewedLanes)) {
        using vecCore::MaskedAssign;
#ifdef DRIVER_PRINT_PROGRESS
        if (partDebug) {
          cout << " --Inserting New Track - part 2/2: (New) 'masked' reset of remaining state: " << endl;
          cout << " *** Existing values - values before change" << endl;
          ReportRowOfDoubles("x1", x1);
          ReportRowOfDoubles("x", x);
          ReportRowOfDoubles("hTry", hTry);
          ReportRowOfDoubles("hdid", hdid);
          ReportRowOfDoubles("numStep", nstp);
          cout << " *** Existing values - recall" << endl;
          ReportRowOfBools<Real_v>("renewedLanes", renewedLanes);
          ReportRowOfDoubles("xStart", xStartLane);
        }
#endif

        // -- Reset working variables for lanes -- vectorised & clean
        MaskedAssign( nstp, renewedLanes, Index_v(0)); // ==> Requires compatible Integer type ...
        // MaskedAssign(  nstp,  renewedLanes, Real_v(0.0) );
        MaskedAssign( x1,   renewedLanes, xStartLane);
        MaskedAssign( x,    renewedLanes, xStartLane);
        MaskedAssign( hTry, renewedLanes, hStepLane);
        // MaskedAssign(  x,  renewedLanes, x1         );       // Done at top of loop - for now
        MaskedAssign( hdid, renewedLanes, Real_v(0.0) ); // Maybe not really needed ... but good!)

        x2 = x1 + hStepLane; // It does not hurt to do it again for some lanes ...

        // Copy the remaining values ...
        // for (unsigned int i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
        //    y[i] = vecCore::Blend( renewedLanes, yInput[i], yNext[i] );
        //    vecCore::MaskedAssign( y, !renewedLanes, yNext[i] );
        // }

        // vecCore::MaskedAssign( isLastStepLane, renewedLanes, Bool_v(false) );
        isLastStepLane = isLastStepLane && !renewedLanes;

#ifdef DRIVER_PRINT_PROGRESS
        MaskedAssign( momStart, renewedLanes, GetMomentumMag(y) );
        
        if (partDebug) {
          cout << " --Inserting New Track - After 'masked' reset of remaining state: " << endl;
          cout << " *** Existing values - values after load / change" << endl;
          ReportRowOfDoubles("x1", x1);
          ReportRowOfDoubles("x", x);
          ReportRowOfDoubles("hTry", hTry);
          ReportRowOfDoubles("hdid", hdid);
          ReportRowOfDoubles("numStep", nstp);
          cout << " *** Existing values - recall" << endl;
          ReportRowOfBools<Real_v>("renewedLanes", renewedLanes);
          ReportRowOfDoubles("xStart", xStartLane);
        }
#endif
        
      }
      
    } // end if ( ! vecCore::MaskFull( laneContinues ) )  // At least one lane is finished

    /*    Bool_v leftLanes = (nstp<=fMaxNoSteps) && (x < x2) && (!isLastStepLane) ;
        int countLeftLanes=0;
        int indLastLane;
        // cout << " leftLanes is: " << leftLanes << std::endl;
        if( !vecCore::MaskEmpty(leftLanes) )
        {
          for (int i = 0; i < vecCore::VectorSize<Real_v>(); ++i)
          {
            if (leftLanes[i] == 1)
            {
              countLeftLanes++;
              indLastLane = i;
              // cout << indLastLane << std::endl;
            }
          }
        }

        // cout<< "countLeftLanes is: "<<countLeftLanes << std::endl;

        if (countLeftLanes == 1)
        {
          // double hstepOneLane = hStepLane[indLastLane] - hDone[indLastLane];
          vecgeom::Vector3D<double> Pos, Mom;
          for (int i = 0; i < 3; ++i)
           {
             Pos[i] = y[i][indLastLane];
             Mom[i] = y[i+3][indLastLane];
           }
          ScalarFieldTrack y_input(Pos, Mom);
          ScalarFieldTrack y_output(Pos, Mom);
          // y_input.SetCurveLength( hDone[indLastLane] ) ;
          fpScalarDriver->AccurateAdvance(y_input, hstep[ indexArr[indLastLane] ] - hDone[indLastLane], epsilon,
       y_output );

          isDoneLane[indLastLane] == true;
          // Store Output
          double y_output_arr[12];
          y_output.DumpToArray(y_output_arr);
          yOutput[indexArr[indLastLane]].LoadFromArray(y_output_arr);
        }*/

  } // end of while loop

  if( numBadSize > 0) {
     cout << "Copying yIn => yOut for " << numBadSize << " 'bad' tracks - i.e. hStep <= 0.0 " << endl;
     for (int i = 0; i < nTracks; ++i) {
        if( badStepSize[i] ) {
           yOutput[i] = yInput[i];
           // cout << " Fixed  yOutput = yInput  for track " << i << " with hStep= " << hstep[i] << endl;
        }
     }
  }
  
  for (int i = 0; i < nTracks; ++i) {
     FieldTrack ft= yOutput[i];
     Vector3D<double> mom ( ft[3], ft[4], ft[5] );
     if( mom.Mag2() == 0.0 ) {
        cout << " Zero output momentum for track " << i << " with hStep = " << hstep[i] << endl;
     }
  }
  
  return;
} // end of AccurateAdvance (flexible / vector )  .....................

// --------------------------------------------------------------------------
#ifdef QUICK_ADV_ARRAY_IN_AND_OUT
template <class Real_v, class T_Stepper, unsigned int Nvar>
typename Mask<Real_v> SimpleIntegrationDriver<Real_v, T_Stepper, Nvar>::QuickAdvance(Real_v yarrin[], // In
                                                                                     const Real_v dydx[],
                                                                                     Real_v hstep, // In
                                                                                     Real_v yarrout[],
                                                                                     Real_v &dchord_step,
                                                                                     Real_v &dyerr) // In length
{
  std::cerr << "ERROR in SimpleIntegrationDriver::QuickAdvance()" << std::endl;
  std::cerr << "      Method is not yet implemented." << std::endl;

  //            FatalException, "Not yet implemented.");
  dyerr = dchord_step = hstep * yarrin[0] * dydx[0];
  yarrout[0]          = yarrin[0];
  exit(1);
}
#endif

template <class T_Stepper, unsigned int Nvar>
void SimpleIntegrationDriver<T_Stepper, Nvar>::ReportInvalidStepInputs(double hStepArr[], int nTracks)
{
  for (int i = 0; i < nTracks; ++i) {
    double hStep = hStepArr[i];
    if (hStep <= 0.0) {
      if (hStep == 0.0) {
        std::cerr << "Proposed step of track " << i << " is zero; hstep = " << hStep << " !";
      } else {
        std::cerr << "Invalid run condition." << std::endl
                  << "Proposed step is negative; (i= " << i << " ) hstep = " << hStep << "." << std::endl;
      }
    }
  }
}

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
bool SimpleIntegrationDriver<T_Stepper, Nvar>::CheckOutput(Real_v Output[], int lane, int initialSlot,
                                                           std::string testName, std::string varName)
{
  bool allGood = true;
  for (unsigned int j = 0; j < Nvar; j++) {
    double current  = vecCore::Get(Output[j], lane);
    double expected = 10. * initialSlot + j;
    double diff     = current - expected;
    if (std::fabs(diff) > 1.e-9 * std::fabs(expected)) {
      std::cerr << testName << " : ERROR in Output " << varName << " [lane= " << lane << " ] "
                << " [iVar= " << j << " ] "
                << " current = " << current << " VS expected = " << expected << "  diff = " << current - expected
                << std::endl;
      allGood = false;
    }
  }
  return allGood;
}

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
bool SimpleIntegrationDriver<T_Stepper, Nvar>::TestInitializeLanes() // int numTracks)
{
  using std::cerr;
  using std::cout;
  using std::endl;
  using vecCore::Get;

  bool allOk = true;

  cout << " TestInitializeLanes() called. " << endl;

  constexpr int numTracks = 8;
  FieldTrack yInput[numTracks];
  bool       badStepSize[numTracks];
  // double xStart[numTracks] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
  double chargeArr[numTracks] = {1.0, -2.0, 3.0, -4.0, 5.1, -6.1, 7.1, -8.1};
  // ------------------------ Test case 1 ------------------------
  Real_v yWorkLanes[Nvar], hStepLane, chargeLane, sBegin;
  double hStep[numTracks] = {0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1};
  //     ======
  int numBadSize = 0;

  cout << "    calling CreateInput() " << endl;
  CreateInput(yInput, numTracks);

  cout << "    test:  chargeArr [0] = " << chargeArr[0] << " ,  [1] = " << chargeArr[1] << ",  [2] = " << chargeArr[2]
       << ",  [3] = " << chargeArr[3] << endl;
  cout << "    calling InitializeLanes() " << endl;
  int nFilled1 = -1;
  int numNext1 = InitializeLanes(yInput, hStep, chargeArr, /* xStart, */ numTracks, 
                                 yWorkLanes, hStepLane, chargeLane, sBegin, nFilled1,
                                 badStepSize, numBadSize);
  cout << " - Init Lanes returned = " << nFilled1 << endl;

  Real_v xLanes(0.0);
  // cout << " - Printing  x, charge, 0.0,  y[] " << endl;
  ReportStatus(xLanes, chargeLane, hStepLane, Real_v(0.0), yWorkLanes);

  std::string testName("Test of InitializeLanes # ");
  int testId = 1;

  const int minSize = std::min((int)vecCore::VectorSize<Real_v>(), (int)numTracks);
  if (nFilled1 != minSize) { // Can cope with 2, 4 or 8
    cerr << testName << testId << " ERROR> Found less than all lanes filled: Number = " << nFilled1
         << " != expected = " << vecCore::VectorSize<Real_v>() << endl;
  }
  assert(nFilled1 == vecCore::VectorSize<Real_v>());

  if (numNext1 != vecCore::VectorSize<Real_v>()) { // Can cope with 2, 4 or 8
    cerr << testName << testId << " ERROR> Found next track number = " << numNext1
         << " != expected = " << vecCore::VectorSize<Real_v>() << endl;
  }
  assert(numNext1 == vecCore::VectorSize<Real_v>());

  if (numBadSize != 0) {
    cerr << testName << testId << " ERROR> Found non-zero bad size lanes." << endl;
  }
  assert(numBadSize == 0);

  // Calling Checking method
  cout << "      calling Checking method" << endl;
  std::string varNameOutput("yWorkLanes");
  for (unsigned int iSlot = 0; iSlot < vecCore::VectorSize<Real_v>(); iSlot++) {
    bool ok = CheckOutput(yWorkLanes, iSlot, iSlot, testName, varNameOutput);
    allOk   = allOk && ok;
  }
  cout << " END of first test of Initialize Lanes. " << endl;
  cout << endl;

  // ------------------------ Test case 2 ------------------------
  Real_v yWorkLanes2[Nvar]; // hStepLane2, sBegin2; ==> To check clean (not overwrite)
  testId = 2;
  cout << " START of 2nd test of Initialize Lanes. " << endl;

  // int identityNum[9]=     { 0,     1,   2,    3,   4,    5,   6,    7,  8 };
  double hStep2[numTracks] = {0.0, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0};
  //     ======
  int nFilled2      = -1;
  int nextLocation2 = InitializeLanes(yInput, hStep2, chargeArr, /* xStart, */ numTracks, 
                                      yWorkLanes2, hStepLane, chargeLane, sBegin, nFilled2,
                                      badStepSize, numBadSize);

  ReportStatus(xLanes, chargeLane, hStepLane, Real_v(0.0), yWorkLanes);

  // int identityNum[9]= { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
  int predictedLast2[9] = {1, 3, 5, 7, 8, 8, 8, 8, 8};
  int expectNext2       = predictedLast2[std::min((int)vecCore::VectorSize<Real_v>(), 4)];
  if (nextLocation2 != expectNext2) {
    cerr << testName << testId << " ERROR> Found less than the expect id of next location: Number = " << nFilled2
         << "  expected = " << expectNext2 << endl;
  }
  int expectFill2 = std::min((int)vecCore::VectorSize<Real_v>(), 3);
  if (nFilled2 != expectFill2) {
    cerr << testName << testId << " ERROR> Found less than the expect number of lanes filled: Number = " << nFilled2
         << "  expected = " << expectFill2 << endl;
  }
  // assert( nFilled2 == std::min( vecCore::VectorSize<Real_v>(), 3 ) );
  // int identityNums2[9]= { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
  int predictedNumBad2[9] = {0, 2, 3, 5, 5, 5, 5, 5, 5};
  int expectedNumBad2 =
      predictedNumBad2[vecCore::VectorSize<Real_v>()]; // ==? vecCore::VectorSize<Real_v>() - nFilled2;
  if (numBadSize != expectedNumBad2)
    cerr << testName << testId << " ERROR> Found " << numBadSize << " bad step-size"
         << " versus " << expectedNumBad2 << " expected for VecSize = " << vecCore::VectorSize<Real_v>() << endl;

  int lane      = 0;
  double hLane0 = vecCore::Get(hStepLane, lane);
  if (hLane0 != 2.0) {
    cerr << testName << testId << " ERROR> hStep[ " << lane << " initial lanes." << endl;
  }

  std::string nameHstep("hStep");
  CheckLaneValue(hStepLane, 0, 2.0, testName, nameHstep);
  CheckLaneValue(hStepLane, 1, 4.0, testName, nameHstep);
  CheckLaneValue(hStepLane, 2, 6.0, testName, nameHstep);

  assert(std::fabs(Get(hStepLane, lane = 0) - 2.0) < 1.0e-8);
  // double hLane1 = Get( hStepLane, 1 );
  assert(std::fabs(Get(hStepLane, lane = 1) - 4.0) < 1.0e-8);
  if (vecCore::VectorSize<Real_v>() > 2) {
    assert(std::fabs(Get(hStepLane, lane = 2) - 6.0) < 1.0e-8);
  }

  int slot = -1;
  bool ok1, ok2, ok3;
  ok1   = CheckOutput(yWorkLanes2, lane = 0, slot = 2, testName, varNameOutput);
  ok2   = CheckOutput(yWorkLanes2, lane = 1, slot = 4, testName, varNameOutput);
  ok3   = CheckOutput(yWorkLanes2, lane = 2, slot = 6, testName, varNameOutput);
  allOk = allOk && ok1 && ok2 && ok3;

  return allOk;
}


template <class T_Stepper, unsigned int Nvar>
inline void SimpleIntegrationDriver<T_Stepper, Nvar>
    ::ComputeAndSetErrcon()
{
  fErrcon = Math::Pow(fMaxSteppingIncrease / fSafetyFactor, 1.0 / GetPowerGrow() );

  std::cout << "SimpleIntegrationDriverComputAndSetErrcon():  fErrcon = " << fErrcon
            << "  from:  maxStepIncrease =  " << fMaxSteppingIncrease
            << "  fSafetyFactor = " << fSafetyFactor
            << "  power-grow =  " << GetPowerGrow() << std::endl;
  
  assert( fErrcon > 0.0 ); 
  // return fErrcon;
}

template <class T_Stepper, unsigned int Nvar>
void SimpleIntegrationDriver<T_Stepper, Nvar>::AccurateAdvance(const FieldTrack yInput[],
                                                               const double     hstep[],
                                                               const double     charge[],
                                                               // double           epsilon, // Can be scalar or varying
                                                               FieldTrack       yOutput[],
                                                               int              nTracks,
                                                               bool             stillOK[] ) const
{
  AccurateAdvance<geant::Double_v>(yInput, hstep, charge,
                                   // epsilon, // Can be scalar or varying
                                   yOutput, stillOK, nTracks);
}

#ifdef EXTEND_SINGLE
template <class T_Stepper, unsigned int Nvar>
void SimpleIntegrationDriver<T_Stepper, Nvar>::AccurateAdvance(const FieldTrack &yInput, const double hstep,
                                                               const double charge,
                                                               // double epsilon,
                                                               FieldTrack &yOutput,
                                                               bool succeeded) const
{
  AccurateAdvance<double>(&yInput, &hstep, &charge,
                          // epsilon, // Can be scalar or varying
                          &yOutput, &succeeded, 1);
}
#endif


#endif /* SimpleIntegrationDriver_Def */
