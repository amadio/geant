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

#ifndef RollingIntegrationDriver_Def
#define RollingIntegrationDriver_Def

// #include "Geant/TemplateFieldTrack.h"
#include "base/AlignedBase.h"
#include "Geant/FieldTrack.h"

// #include "TemplateVScalarIntegrationStepper.h"
// #include "IntegrationStepper.h"

#include "base/Vector.h"

// Adding because adding scalar stepper for new constructor (KeepStepping)
// #include "Geant/VScalarIntegrationStepper.h"

// Adding to send in scalar driver to deal with 1/2 remaining lanes
// #include "Geant/ScalarIntegrationDriver.h"
// #include "Geant/ScalarFieldTrack.h"

#include "Geant/FlexIntegrationDriver.h"

#include "Geant/BaseRkIntegrationDriver.h"
#include "Geant/AuxVecMethods.h"

#define  USE_ERROR_ESTIMATOR 1

#ifdef USE_ERROR_ESTIMATOR
#include "ErrorEstimatorSixVec.h"
#endif
    
// #include "Geant/FormattedReporter.h"  // Direct include is not needed.
#include "Geant/PrintDriverProgress.h"

#include "Geant/VectorTypes.h" //  Defines geant::Double_v
#include "Geant/math_wrappers.h"

#ifndef NO_FIELD_STATISTICS
#define GVFLD_STATS 1
#endif

// --------------------------------------------------------------
template <class T_Stepper, unsigned int Nvar>
class RollingIntegrationDriver :
   // public FlexIntegrationDriver,
   public BaseRkIntegrationDriver<T_Stepper, Nvar>,
   public vecgeom::AlignedBase {
public:
  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

  RollingIntegrationDriver(double     hminimum, // same
                          T_Stepper *pStepper,
                          int        numberOfComponents = 6,
                          int        statsVerbosity = 1);

  ~RollingIntegrationDriver();

  virtual void AccurateAdvance(const FieldTrack yInput[], const double hstep[], const double charge[], double epsilon,
                               FieldTrack yOutput[], int nTracks, bool succeeded[]) const override final;

#ifdef EXTEND_SINGLE
  virtual void AccurateAdvance(const FieldTrack &yInput, const double hstep, const double charge, double epsilon,
                               FieldTrack &yOutput, bool succeeded) const override final;
#endif

  // Implemented in terms of the following templated function:
  template <class Real_v>
  void AccurateAdvance(const FieldTrack yInput[], const double hstep[], const double charge[], double epsilon,
                       FieldTrack yOutput[], bool succeeded[], int nTracks) const;
  // Drive Runge-Kutta integration of ODE for several tracks (ntracks)
  // with starting values yInput, from current 's'=0 to s=h with variable
  // stepsize to control error, so that it is bounded by the relative
  // accuracy eps.  On output yOutput is value at end of interval.
  // The concept is similar to the odeint routine from NRC 2nd edition p.721

  // unsigned long GetNumberOfStepperCalls() { return fStepperCalls; }
  // unsigned long GetNumberOfTotalSteps()   { return fNoTotalSteps; }

  /*****
      inline void   GetDerivatives( const TemplateFieldTrack<Real_v> &y_curr,     // const, INput
                                          Real_v    charge,
                                          Real_v    dydx[]   );  //       OUTput
   ******/

  // EquationOfMotion<Real_v>* GetEquationOfMotion() { return fpStepper->GetEquationOfMotion(); }
  // const EquationOfMotion<Real_v>* GetEquationOfMotion() const { return fpStepper->GetEquationOfMotion(); }

  // RollingIntegrationDriver* Clone() const;
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

/******
  template <class Real_v>
  Real_v ComputeNewStepSize(Real_v errMaxNorm,    // normalised error
                            Real_v hstepCurrent); // current step size
                                                  // Taking the last step's normalised error, calculate
                                                  // a step size for the next step.
                                                  // Do not limit the next step's size within a factor of the
                                                  // current one.

  template <class Real_v>
  Real_v ComputeNewStepSize_WithinLimits(Real_v errMaxNorm,    // normalised error
                                         Real_v hstepCurrent); // current step size
                                                               // Taking the last step's normalised error, calculate
                                                               // a step size for the next step.
  // Limit the next step's size within a range around the current one.
*****/
  
  template <class Real_v>
  int InitializeLanes(const FieldTrack yInput[], const double hstep[], const double charge[], int nTracks,
                      //         bool        badStepSize[],   // Output - redundant
                      int indexArr[], // [vecCore::VectorSize<Real_v>()] - Output
                      Real_v y[],     // [Nvar]        - Output
                      Real_v &hStepLane, Real_v &chargeLane, Real_v &startCurveLength,
                      int &numFilled, // [Out]: number of filled lanes
                      int &numBadSize) const;
  // Load 'work' into array y, indices of lanes into indexArr, etc
  // Return array index of 'next' track - i.e. first track not yet loaded

  template <class Real_v>
  bool InsertNewTrack(const FieldTrack yInput[], const double hstep[], const double charge[], const int slot,
                      int &trackNextInput, bool succeeded[], Real_v y[], Real_v &hStepLane, Real_v &chargeLane,
                      Real_v &startCurveLength,
                      int indexArr[], // [vecCore::VectorSize<Real_v>()]
                      int nTracks) const;
  template <class Real_v>
  void StoreOutput(const Real_v yEnd[], const Real_v x,
                   int currIndex, // Index in Real_v
                   FieldTrack yOutput[],
                   int indOut, // location in Ouptut
                   const double hstep[], bool succeeded[], int nTracks) const;

  void ComputeAndSetErrcon() { BaseRkIntegrationDriver<T_Stepper, Nvar>::ComputeAndSetErrcon(); }
  void CheckParameters()     { BaseRkIntegrationDriver<T_Stepper, Nvar>::CheckParameters(); }

  // Access parameters
  double GetPowerShrink()  const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fPowerShrink; }
  double GetPowerGrow()    const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fPowerGrow; }
  double GetSafetyFactor() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fSafetyFactor; }
  int    GetVerboseLevel() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fVerboseLevel; }
  int    GetStatisticsVerboseLevel() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStatisticsVerboseLevel(); }
  unsigned int GetStepperOrder() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStepperOrder(); }

  
  double GetMinimumStep() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetMinimumStep(); }
  double GetErrcon()      const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetErrcon(); }

  // Modifier methods
  int IncrementStepperCalls() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::IncrementStepperCalls() ; }
  void SetMaxNoSteps(int val) { BaseRkIntegrationDriver<T_Stepper,Nvar>::SetMaxNoSteps(val) ; }
  
  /**
  const T_Stepper *GetStepper() const { return fpStepper; } // { BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStepper(); }
        T_Stepper *GetStepper()       { return fpStepper; }   // { BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStepper(); }
   **/
  void CreateInput(FieldTrack yInput[], int nTracks)
  { BaseRkIntegrationDriver<T_Stepper,Nvar>::CreateInput( yInput, nTracks);  }
        
// BaseRkIntegrationDriver<T_Stepper,Nvar>::  
   
public: // For now
  template <class Real_v>
  bool TestInitializeLanes(); // (int numTracks)
                              // Simple, unit-test like check.  Returns ok or not
  template <class Real_v>
  bool CheckOutput(Real_v yOutput[], int lane, int initialSlot, std::string testName, std::string varName);
  // returns 'true' if ok

  static constexpr int fMaxStepBase = 250;  
  static constexpr double fMaxSteppingIncrease = 5.0;
  static constexpr double fMaxSteppingDecrease = 0.1;
  
private:
  // Private copy constructor and assignment operator.

  RollingIntegrationDriver(const RollingIntegrationDriver &);
  // Copy constructor used to create Clone method

  RollingIntegrationDriver &operator=(const RollingIntegrationDriver &) = delete;

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

  // double fMinimumStep; // same
  // Minimum Step allowed in a Step (in absolute units)
  static constexpr double fSmallestFraction = 1.0e-7; // Expected value: larger than 1e-12 to 5e-15;

  const double fHalfPowerShrink;
  unsigned long fMaxNoSteps= 100;  // Default - expect it to be overriden in constructor
  
  // ---------------------------------------------------------------
  // Compilation constants
  const bool partDebug  = false;                 // Enforce debugging output
  const bool useOneStep = true;                  //  Algorithm selection - false for KeepStepping

  // ---------------------------------------------------------------
  //  STATE
  // public:

  // Related to AccurateAdvance
  // Variables required for track insertion algorithm
  // int    fNTracks = 0;        //  number under current integration -> ~ current array size
  // ScalarIntegrationStepper *fpScalarStepper= nullptr;  // --> call Stepper with scalar values (args)
  // ScalarIntegrationDriver  *fpScalarDriver = nullptr;  // --> (later) use this driver with scalar args
  // bool partDebug = false ;
  // ---

  mutable unsigned long fNoTotalSteps = 0; // , fNoBadSteps = 0, fNoSmallSteps = 0, fNoInitialSmallSteps = 0;

}; // End of class definition -- RollingIntegrationDriver

#ifdef DRIVER_PRINT_PROGRESS  
#include "Geant/FormattedReporter.h"
#endif

/*****
template <class Real_v, class T_Stepper, unsigned int Nvar>
constexpr double RollingIntegrationDriver<Real_v, T_Stepper, Nvar>::fMaxSteppingIncrease;

template <class Real_v, class T_Stepper, unsigned int Nvar>
constexpr double RollingIntegrationDriver<Real_v, T_Stepper, Nvar>::fMaxSteppingDecrease;
 *****/


/*********
// ---------------------------------------------------------

template <class Real_v, class T_Stepper, unsigned int Nvar>
inline
void RollingIntegrationDriver<Real_v, T_Stepper, Nvar>
  ::GetDerivatives(const TemplateFieldTrack<Real_v> &y_curr, // const, INput
                         Real_v  charge,
                         Real_v  dydx[])  // OUTput
{
  Real_v  tmpValArr[ncompSVEC]; // TemplateFieldTrack<Real_v>::ncompSVEC];
  y_curr.DumpToArray( tmpValArr  );
  fpStepper -> RightHandSideVIS( tmpValArr , charge, dydx );
}
 ********/

// template <class T_Stepper, unsigned int Nvar>
// const int  RollingIntegrationDriver< /* Real_v, */ T_Stepper, Nvar>::fMaxStepBase = 250;  // Was 5000

// To add much printing for debugging purposes, uncomment the following
// and set verbose level to 1 or higher value !
// #define  GUDEBUG_FIELD 1

// ---------------------------------------------------------

//  Constructor
//
template <class T_Stepper, unsigned int Nvar>
RollingIntegrationDriver<T_Stepper, Nvar>::RollingIntegrationDriver(double     hminimum,
                                                                  T_Stepper *pStepper,
                                                                  int        numComponents,
                                                                  int        statisticsVerbose)
    :
      BaseRkIntegrationDriver<T_Stepper, Nvar>(hminimum,
                                               pStepper,     // or pass pStepper->GetIntegratorOrder()  ??
                                               numComponents,
                                               statisticsVerbose),
      fHalfPowerShrink( 0.5 * BaseRkIntegrationDriver<T_Stepper, Nvar>::GetPowerShrink() )
      // , 
      // fMinimumStep(hminimum),
      // -- fSmallestFraction( 1.0e-12 ),
      // -- fNoIntegrationVariables(numComponents),  // ==> Nvar
      // -- fMinNoVars(12), fNoVars(std::max((int)Nvar, std::max((int)fMinNoVars, (int)numComponents))),
      // - fPowerShrink(-1.0 / pStepper->GetIntegratorOrder()),       //  exponent for shrinking
      // -fPowerGrow(-1.0 / (1.0 + pStepper->GetIntegratorOrder())), //  exponent for growth
      // - fErrcon(0.0),
      // fStatisticsVerboseLevel(statisticsVerbose),
      // -- fNoInitialSmallSteps(0)
      // fVerboseLevel(0)
{
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
  assert(pStepper != nullptr);
  assert(Nvar <= (unsigned int)numComponents /* was ncompSVEC */ ); // Ensure that arrays are large enough for Integr.

  // fpStepper = pStepper;

  // BaseRkIntegrationDriver<T_Stepper, Nvar>::ComputeAndSetErrcon();
  ComputeAndSetErrcon();
  SetMaxNoSteps( fMaxStepBase / pStepper->GetIntegratorOrder() );

  BaseRkIntegrationDriver<T_Stepper, Nvar>::ComputeAndSetErrcon();  
  // ComputeAndSetErrcon();

  CheckParameters();

#ifdef GUDEBUG_FIELD
  fVerboseLevel = 2;
#endif

  if (GetVerboseLevel()) {
    std::cout << "SiD:ctor> Stepper Order= " << pStepper->GetIntegratorOrder() << " > Powers used: "
              << " shrink = " << GetPowerShrink() << "  grow = " << GetPowerGrow() << std::endl;
  }
  if ((GetVerboseLevel() > 0) || (GetStatisticsVerboseLevel() > 1)) {
     std::cout << "RollingIntegrationDriver created. " << std::endl;
     //         << "invE_nS, QuickAdv-2sqrt with Statistics " << fStatsStatus << std::endl;
     // ( fStatsEnabled ? "enabled" : "disabled" )
  }

  // For track insertion
}

// ---------------------------------------------------------

//  Copy Constructor - used by Clone
template <class T_Stepper, unsigned int Nvar>
RollingIntegrationDriver<T_Stepper, Nvar>::RollingIntegrationDriver(
    const RollingIntegrationDriver</*Real_v,*/ T_Stepper, Nvar> &right)
    :
   BaseRkIntegrationDriver<T_Stepper, Nvar>::BaseRkIntegrationDriver( right )

// fMinimumStep(right.fMinimumStep), fSmallestFraction(right.fSmallestFraction),
      // fNoIntegrationVariables( right.fNoIntegrationVariables ),
      // fMinNoVars(right.fMinNoVars), fNoVars(std::max((int)Nvar, fMinNoVars)), 
      // fPowerShrink(right.fPowerShrink), fPowerGrow(right.fPowerGrow), fErrcon(right.fErrcon),
      // fStatisticsVerboseLevel(right.fStatisticsVerboseLevel) // ,
      // fVerboseLevel(right.fVerboseLevel)
{
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
   
  // const T_Stepper *protStepper = right.fpStepper;
  // fpStepper                    = protStepper->Clone();
  // ==> This should happen in the base class .....
   
  ComputeAndSetErrcon();
  fMaxNoSteps = std::max( fMaxStepBase / GetStepperOrder(), 10 );  // fpStepper->GetIntegratorOrder();

  if ((GetVerboseLevel() > 0) || (GetStatisticsVerboseLevel() > 1)) {
    std::cout << "RollingIntegrationDriver copy Constructor. "
       // << "invE_nS, QuickAdv-2sqrt with Statistics "
       // << fStatsStatus
              << std::endl;
  }
}

// ---------------------------------------------------------

//  Destructor
template <class T_Stepper, unsigned int Nvar>
RollingIntegrationDriver<T_Stepper, Nvar>::~RollingIntegrationDriver()
{
   if (GetStatisticsVerboseLevel() > 1) {
      // PrintStatisticsReport();
   }

   // delete fpScalarDriver;
   // delete fpScalarStepper;
   // delete fpStepper;
}

// ---------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void RollingIntegrationDriver<T_Stepper, Nvar>::OneGoodStep(const Real_v yStart[], const Real_v dydx[],
                                                           const Real_v charge,
                                                           Real_v &x, // InOut
                                                           Real_v htry,
                                                           double eps_rel_max,
                                                           // const Real_v  eps_rel_max,
                                                           Real_v yFinal[], // Out-values
                                                           Real_v &hdid,    // Out - achieved length
                                                           Real_v &hnext)   // Out - proposed next integration length
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

  using vecCore::math::Min;
  using vecCore::math::Max;
  using vecCore::math::Exp;
  using vecCore::math::Log;
  using vecCore::Get;

  using FormattedReporter::ReportRowOfDoubles;
  using FormattedReporter::ReportRowOfSquareRoots;
  using FormattedReporter::ReportManyRowsOfDoubles;
  using FormattedReporter::ReportRowOfBools;
  using ReportValuesOfVectors::ReportConditionLanes;

  if (partDebug) { cout << "\n" << endl; }
  // const int vecCore::VectorSize<Real_v>() = vecgeom::vecCore::VectorSize<Real_v>();
  // const int ncompSVEC = TemplateFieldTrack<Real_v>::ncompSVEC;

  Real_v errmaxSqFallThru(0.0);
  
  Real_v xnew, yerr[Nvar /*ncompSVEC*/], ytemp[Nvar /*ncompSVEC*/];
  Real_v h = htry; // Set stepsize to the initial trial value
  // Renamed it to hStep

  // Real_v  //  Epsilon was variable per track (ToCheck)
  // double invEpsilonRelSq = 1.0 / (eps_rel_max * eps_rel_max);

  static int tot_no_trials = 0; // Should be thread_local - or suppressed. Just statistics
  const int max_trials     = 100;

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

  do {
    // Bool_v alreadyFinished = finished;  // State at start of iteration
    Bool_v Active = !finished;
    Real_v errmax_sq= -1.0;

    itersLeft--;
    iter++;

    if (partDebug) std::cout << " OneGoodStep - iteration = " << iter << endl;

    // #ifdef STORE_ONCE
    vecCore::MaskedAssign(h, finished, Real_v(0.0)); // Set h = 0.0 for finished lanes -- ensure no change !
                                                     // #endif

    // if ( !vecCore::IsFull(stepSizeUnderflow || goodStep) )
    // {
    fpStepper->StepWithErrorEstimate(yStart, dydx, charge, h, ytemp, yerr); // CAREFUL -> changes for others ?
   
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

#ifdef USE_ERROR_ESTIMATOR
    ErrorEstimatorSixVec fErrorEstimator( eps_rel_max, GetMinimumStep() );
    errmax_sq = fErrorEstimator.EstimateError( yerr, h, ytemp );
#else
    Real_v errpos_sq = 0.0; // square of displacement error
    Real_v errmom_sq = 0.0; // square of momentum vector difference
    Real_v epsPosition = eps_rel_max * vecCore::math::Max(h, Real_v(fMinimumStep)); // Uses remaining step 'h'
    // Could change it to use full step size ==> move it outside loop !!   2017.11.10 JA
    Real_v invEpsPositionSq = 1.0 / (epsPosition * epsPosition);

    // Evaluate accuracy
    Real_v errpos_sq = yerr[0] * yerr[0] + yerr[1] * yerr[1] + yerr[2] * yerr[2];
    errpos_sq * invEpsPositionSq; // Scale relative to required tolerance

    // Accuracy for momentum
    Real_v magmom_sq = ytemp[3] * ytemp[3] + ytemp[4] * ytemp[4] + ytemp[5] * ytemp[5];
    Real_v sumerr_sq = yerr[3] * yerr[3] + yerr[4] * yerr[4] + yerr[5] * yerr[5];

    // vecCore::CondAssign(magmom_sq > 0.0, sumerr_sq/magmom_sq, sumerr_sq, &errmom_sq);
    constexpr double tinyValue = 1.0e-80; // Just to ensure there is no division by zero
    Real_v errmom_sq           = sumerr_sq / (magmom_sq + tinyValue);

    errmom_sq *= invEpsilonRelSq;
    errmax_sq = vecCore::math::Max(errpos_sq, errmom_sq); // Square of maximum error
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

    Bool_v laneDone = (goodStep | finished);
    bool allDone    = vecCore::MaskFull(laneDone);

    finished = laneDone;
    Active   = !finished;

#ifdef DRIVER_PRINT_PROGRESS
    if (partDebug) {
      ReportRowOfBools<Real_v>("goodStep", goodStep);
      ReportRowOfBools<Real_v>("laneDone", laneDone);
      ReportRowOfBools<Real_v>("(updated) finished", finished);
    }
#endif
    
    if (allDone) // All (or remaining) steps succeeded.
    {
      // Idea 1.5
      if (partDebug) cout << "Store and Report Stored lanes - v1.5 allDone - about to break." << endl;

      StoreGoodValues(ytemp, h, errmax_sq, goodStep, yFinal, hFinal, errmax_sqFinal);

      break;
    }

    Real_v errPower = PowerIf<Real_v>(errmax_sq,
                                      fHalfPowerShrink, // 0.5 * GetPowerShrink(),
                                      !laneDone);
    Real_v hReduced = h * Max(Real_v(0.1), GetSafetyFactor() * errPower);

    Real_v hnew = vecCore::Blend(finished, Real_v(0.0), hReduced);
    xnew        = x + hnew;

    stepSizeUnderflow = Active && (xnew == x);

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
      StoreGoodValues(ytemp, h, errmax_sq,
                      (goodStep || stepSizeUnderflow), // && !alreadyFinished,
                      yFinal, hFinal, errmax_sqFinal);
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

    Real_v errmaxSqFallThru(0.0);
    Real_v errmaxSqFallThru= errmax_sq;    
    
  } while (itersLeft > 0 && (!vecCore::MaskFull(finished)) //  was MaskFull( stepSizeUnderflow || goodStep ) )
           );

  tot_no_trials += iter;

#ifdef STORE_ONCE
  //  'Idea 3' - Store exactly one time ( here - except on loop exit)
  StoreGoodValues(ytemp, h, errmaxSqFallThru, finished, yFinal, hFinal, errmax_sqFinal);
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

  // #ifdef CHECK_STRETCH_FACTOR
  // ------------------------------------------
  // The old way (improved) - to cross check
  constexpr double minErr2 = 1e-100;
  Real_v emax2pos          = Max(errmax_sqFinal, Real_v(minErr2));
  Real_v errStretchOld     = GetSafetyFactor() * Exp((0.5 * GetPowerGrow()) * Log(emax2pos)); // Was: Log(errmax_sqFinal) );
  // ReportRowOfDoubles( "-raw-errStretch", errStretch);
  errStretchOld  = Min(errStretchOld, Real_v(fMaxSteppingIncrease));
  Bool_v zeroErr = errmax_sq <= minErr2;
  // Fix the size for steps with zero error !!
  vecCore::MaskedAssign(errStretchOld, zeroErr, Real_v(fMaxSteppingIncrease));
  // ReportRowOfDoubles( "old: errStretch", errStretchOld);
  // ------------------- End of Old way -----
  // #endif

  // Check against fErrcon to avoid calling power ... saves work if any are 'over' max
  auto errcon = GetErrcon();
  Bool_v underThresh = errmax_sqFinal <= errcon * errcon;
  Real_v errStretch1raw  = fSafetyFactor * PowerIf(errmax_sqFinal, 0.5 * fPowerGrow, !underThresh);
  // Note:  lanes with 'false' argument (i.e. underThresh=true) will have value 1.0
  Real_v errStretch =
    vecCore::Blend( underThresh, Real_v(fMaxSteppingIncrease), errStretch1raw );

  hnext = errStretch * h;

#ifdef CHECK_STRETCH_FACTOR
  if (!vecCore::MaskEmpty(errStretch - errStretchOld > 1.0e-12 * errStretch)) {
    cout << "ERROR> Lanes found with differences in calculated value of 'errStretch'"
         << "       ( factor for stretching step size for 'good' next step." << endl;
    ReportRowOfDoubles("old-new: errStretch", errStretch - errStretchOld);
    ReportRowOfDoubles("old: errStretch", errStretchOld);
    ReportRowOfDoubles("new: errStretch", errStretch);
  }
#endif

  // ReportRowOfDoubles( "OGS: h-final", hFinal);

  hdid = hFinal;
  x += hdid;

#ifdef DRIVER_PRINT_PROGRESS  
  bool OGSreport = true;
  if (partDebug && OGSreport) {
    ReportRowOfDoubles("OGS: errmax2", errmax_sq);
    ReportRowOfDoubles("OGS: h-did ", hdid);
    ReportRowOfDoubles("OGS:  new-x", x);
    ReportRowOfDoubles("OGS: h-next", hnext);
    ReportRowOfDoubles("OGS: facStretch", errStretch);
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
inline void RollingIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::StoreGoodValues(
    const Real_v yWork[], const Real_v &hValue, const Real_v &errMaxSqValue, const vecCore::Mask_v<Real_v> &storeFlag,
    Real_v yFinal[], Real_v &hFinalStore, Real_v &errMaxSqStore) const
// yWork,  hValue,      epsSqVal   represent the output variables
// yFinal, hFinalStore, epsSqStore represent the output variables
{
  using std::cout;
  using std::endl;
  using vecCore::MaskedAssign;
  using vecCore::Get;
  using vecCore::Set;

  if (vecCore::MaskFull(storeFlag)) {
    // std::cout << "StoreGoodValues: Unconditional assignment to output - all together." << std::endl;

    for (unsigned int j = 0; j < Nvar; ++j)
      yFinal[j]         = yWork[j];

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
    using FormattedReporter::ReportRowOfBools;
    using FormattedReporter::ReportRowOfDoubles;
    using FormattedReporter::ReportManyRowsOfDoubles;
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
int RollingIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::InitializeLanes(
    const FieldTrack yInput[], const double hstep[], const double charge[],
    // const  double     xStart [], // --> Now in FieldTrack
    int nTracks,
    //     bool       badStepSize[],  // Output - redundant
    int indexArr[], // [vecCore::VectorSize<Real_v>()]
    Real_v y[], Real_v &hStepLane, Real_v &chargeLane, Real_v &startSlen,
    int &numFilled, // How many were loaded.
    int &numBadSize) const
// Converts input scalar stream to acceptable form of Vc vectors
// for vector processing in OneStep
{
  using vecCore::Set;
  using FormattedReporter::ReportArray;
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

  numBadSize  = 0; // Ensure it starts at Zero !?
  int j       = 0;
  size_t slot = 0;
  do {
    double hVal       = hstep[j];
    double chargeVal  = charge[j];
    bool invalidTrack = (hVal <= 0.0) || (chargeVal == 0.0);
    // badStepSize[j] =  invalidTrack;

    if (invalidTrack) {
      if (hVal <= 0.0)
        std::cout << " WARNING> Non-positive step-size h = " << std::setw(5) << hVal << " for j= " << j << std::endl;
      if (chargeVal == 0.0) std::cout << " WARNING> Zero charge " << chargeVal << " for j= " << j << std::endl;
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
bool RollingIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::InsertNewTrack(
    const FieldTrack yInput[], const double hstep[], const double charge[], const int slot, int &trackNextInput,
    bool succeeded[],
    Real_v y[], // [Nvar]
    Real_v &hStepLane, Real_v &chargeLane, Real_v &startCurveLength,
    int indexArr[], // [vecCore::VectorSize<Real_v>()]
    int nTracks) const

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
    if (hStepNext > 0) {
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

      filled = true;
    } else {
      // A zero or negative step is anomalous - an error
      succeeded[trackNextInput] = (hStepNext == 0.0);
      if (hStepNext == 0) {
        std::cerr << "Proposed step is zero; hstep = " << hStepNext << " !" << std::endl;
      } else {
        std::cerr << "Invalid run condition." << std::endl
                  << "Proposed step is negative; hstep = " << hStepNext << "." << std::endl;
      }
    }
    trackNextInput++;
  }

  return filled;

} // End of InsertNewTrack function

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void RollingIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::StoreOutput(const Real_v yEnd[], const Real_v x,
                                                                       int currIndex, FieldTrack yOutput[], int indOut,
                                                                       const double hstep[], bool succeeded[],
                                                                       int nTracks) const
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
    std::cout << "----Storage position (out-arr): " << indOut
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
      yOutOneArr[i] = vecCore::Get(yEnd[i], currIndex);
    }
    yOutput[indOut].LoadFromArray(yOutOneArr, Nvar);            // fNoIntegrationVariables);
    yOutput[indOut].SetCurveLength(vecCore::Get(x, currIndex)); // x is a double_v variable
  } else {
    succeeded[indOut] = (hOriginal == 0.0);
  }

} // End of StoreOutput function

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void RollingIntegrationDriver<T_Stepper, Nvar>::AccurateAdvance(const FieldTrack yInput[], const double hstep[],
                                                               const double charge[],
                                                               double epsilon, // Can be scalar or varying
                                                               FieldTrack yOutput[], bool stillOK[], int nTracks) const
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

  const std::string methodName = "SID::AccurateAdvance";
#ifdef DRIVER_PRINT_PROGRESS  
  using FormattedReporter::ReportRowOfBools;
  using FormattedReporter::ReportRowOfDoubles;
  using FormattedReporter::ReportRowOfSquareRoots;
  using FormattedReporter::ReportManyRowsOfDoubles;
  using FormattedReporter::ReportRowsOfPositionsMomenta;
  using FormattedReporter::GetMomentumMag;
  using FormattedReporter::ReportArray;
#endif
  using Bool_v = vecCore::Mask_v<Real_v>;
  using vecCore::math::Min;
  using vecCore::math::Max;
  using std::cout;
  using std::endl;

  constexpr unsigned int VecSize = vecCore::VectorSize<Real_v>();
  int indexArr[VecSize]; // vecCore::VectorSize<Real_v>()];

  // Working variables for integration
  Real_v x, hnext, hdid, h, chargeLane, x1, x2, xStartLane, hStepLane;
  Real_v y[Nvar /*ncompSVEC*/ ];     // Working array 1
  Real_v yNext[Nvar /*ncompSVEC*/ ], // Working array 2
      dydx[Nvar /*ncompSVEC*/ ];

  // Real_v ySubStepStart[Nvar /*ncompSVEC*/ ];
  // bool badStepSize[nTracks];      // Step size is illegal, ie. h <= 0

  std::fill_n(stillOK, nTracks, 1); //  Assume success, flag failure ... Strange/JA
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
    indexArr[i]       = -1;

  int idNext = InitializeLanes(yInput, hstep, charge, /* xStart, */ nTracks, /* badStepSize,*/ indexArr, y, hStepLane,
                               chargeLane, xStartLane, numFilled, numBadSize);
  // OLD:
  // int     idNext = vecCore::VectorSize<Real_v>(); // ASSUMES that nTracks >= vecCore::VectorSize<Real_v>() !!! FIX

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

  // Real_v xStart= x;
  h  = hStepLane;
  x1 = xStartLane;
  x2 = x1 + hStepLane;

  x = x1;

  while (
      (   !vecCore::MaskFull(isDoneLane)
          && !vecCore::MaskEmpty((nstp <= fMaxNoSteps)
                                 && (x < x2) && (!isLastStepLane)))
      ||
          idNext < nTracks
     ) {
#ifdef DRIVER_PRINT_PROGRESS     
    if (partDebug)
      std::cout << "************************************" << std::endl
                << "** Top of while loop ***************" << endl
                << "----hStepLane is: " << hStepLane << endl;
#endif
    // if( h > fMinimumStep ) { QuickAdvance .. } else { .. below  //  ( Sequential code  )
    // if (useOneStep) {
    fpStepper->RightHandSideInl(y, chargeLane, dydx); // TODO: change to inline
    //---------****************-----------------------
    OneGoodStep<Real_v>(y, dydx, chargeLane, x, h, epsilon, yNext, hdid, hnext);
    //*********---------------------------------------------------------
    // } else KeepStepping( y, dydx, x, h, epsilon, hdid, hnext, hStepLane, hDone) ;

#ifdef DRIVER_PRINT_PROGRESS    
    if (partDebug) {
      cout << "### Accurate Advance ---- After return from OneGood Step" << endl;
      ReportManyRowsOfDoubles("yStart", y, Nvar);
      ReportManyRowsOfDoubles("dydx", dydx, Nvar);
      ReportRowOfDoubles("h-ask", h);
      ReportRowOfDoubles("h-did", hdid);
      ReportRowOfDoubles("x", x);
      ReportRowOfDoubles("(end) x2", x2);
      cout << "##-------------------------------------------------------------------------------" << endl;
      // ReportManyRowsOfDoubles( "yNext", yNext, Nvar);
      Real_v momStart = GetMomentumMag(y);
      ReportRowsOfPositionsMomenta("yNext", yNext, Nvar, momStart);
    }
#endif

    // lastStepOK = (hdid == h);
    fNoTotalSteps++;

#ifdef DRIVER_PRINT_PROGRESS      
    bool reportMove = true;
    if (partDebug && reportMove) {
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
    Real_v stepThreshold           = vecCore::math::Min(epsilon * hStepLane, fSmallestFraction * (xStartLane + hdid));
    Bool_v avoidNumerousSmallSteps = h < stepThreshold;

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
      cout << " x-Now   = " << x << endl;
      cout << " x2 -x   = " << xremains << endl;
    }
#endif

    hnext = Max(hnext, Real_v(GetMinimumStep()));
    // Note: This has potential for instability i.e. if MinStep is 'too long'

    // Ensure that the next step does not overshoot
    hnext = Min(xremains, hnext);

#ifdef DRIVER_PRINT_PROGRESS      
    if (partDebug) cout << "AccurateAdvance: hnext = " << hnext << " to replace h = " << h << endl;
#endif

    h = hnext;

    // When stepsize overshoots, decrease it!
    // Must cope with difficult rounding-error issues if hstep << x2

    isLastStepLane = (h == 0.0) || isLastStepLane;

    if (partDebug) std::cout << " lastStep : " << isLastStepLane << std::endl;

    nstp += 1; // nstp++;

    succeededLane = (xremains <= 0.0); // (x>=x2); // If it was a "forced" last step ?

    Bool_v laneContinues = (nstp <= fMaxNoSteps) && !succeededLane && !isLastStepLane;

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
        Bool_v CondNoOfSteps = (nstp <= fMaxNoSteps);
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
        if (vecCore::Get(finishedLane, i) && indexArr[i] != -1) {
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
            bool filled = InsertNewTrack(yInput, hstep, charge, i, idNext, stillOK, y, hStepLane, chargeLane,
                                         xStartLane, indexArr, nTracks);
            // isDoneLane[i]   = !filled;
            vecCore::Set(isDoneLane, i, !filled);
            // finishedLane[i] = !filled;
            vecCore::Set(finishedLane, i, !filled);

            Set(renewedLanes, i, filled);

#ifdef DRIVER_PRINT_PROGRESS            
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
          ReportRowOfDoubles("h", h);
          ReportRowOfDoubles("hdid", hdid);
          ReportRowOfDoubles("numStep", nstp);
          cout << " *** Existing values - recall" << endl;
          ReportRowOfBools<Real_v>("renewedLanes", renewedLanes);
          ReportRowOfDoubles("xStart", xStartLane);
        }
#endif
        // -- Reset working variables for lanes -- vectorised & clean
        MaskedAssign(nstp, renewedLanes, Index_v(0)); // ==> Requires compatible Integer type ...
        // MaskedAssign(  nstp,  renewedLanes, Real_v(0.0) );
        MaskedAssign(x1, renewedLanes, xStartLane);
        MaskedAssign(x, renewedLanes, xStartLane);
        MaskedAssign(h, renewedLanes, hStepLane);
        // MaskedAssign(  x,  renewedLanes, x1         );       // Done at top of loop - for now
        MaskedAssign(hdid, renewedLanes, Real_v(0.0)); // Maybe not really needed ... but good!

        x2 = x1 + hStepLane; // It does not hurt to do it again for some lanes ...

        // Copy the remaining values ...
        // for (unsigned int i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
        //    y[i] = vecCore::Blend( renewedLanes, yInput[i], yNext[i] );
        //    vecCore::MaskedAssign( y, !renewedLanes, yNext[i] );
        // }

        // vecCore::MaskedAssign( isLastStepLane, renewedLanes, Bool_v(false) );
        isLastStepLane = isLastStepLane && !renewedLanes;

#ifdef DRIVER_PRINT_PROGRESS        
        if (partDebug) {
          cout << " *** Vectors changed together after 'loop':" << endl;
          ReportRowOfDoubles("x1", x1);
          ReportRowOfDoubles("x", x);
          ReportRowOfDoubles("h", h);
          ReportRowOfDoubles("hdid", hdid);
          ReportRowOfDoubles("numStep", nstp);
          // ReportRowOfDoubles( "hDone",   hDone);
        }
      } else {
        if (partDebug) cout << "-- Insert New Track - part 2/2:  No New Tracks found." << endl;
#endif
      }
      
#ifdef DRIVER_PRINT_PROGRESS        
      // bool ReportAtBottomOfLoop= true;
      // if (ReportAtBottomOfLoop && partDebug)
      bool ReportAfterResets = true;
      if (ReportAfterResets && partDebug) {
        cout << " After all resets --  " << endl;
        cout << " ====================================================================" << endl;
        ReportManyRowsOfDoubles("yCurrent", y, Nvar);
        ReportRowOfDoubles("charge", chargeLane);
        ReportRowOfDoubles("x", x);
        ReportRowOfDoubles("hdid", hdid);
        ReportRowOfDoubles("h(next)", h);
        ReportRowOfBools<Real_v>("isLastStep", isLastStepLane);
        ReportRowOfBools<Real_v>("isDone", isDoneLane);
      }
#endif

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

} // end of AccurateAdvance (flexible / vector )  .....................

//----------------------------------------------------------------------

/*********************************
#define SQR(a)   ((a)*(a))

// QuickAdvance just tries one Step - it does not ensure accuracy
template <class Real_v, class T_Stepper, unsigned int Nvar>//
   typename vecCore::Mask_v<Real_v>
// RollingIntegrationDriver<Real_v, T_Stepper, Nvar>
RollingIntegrationDriver< T_Stepper, Nvar>
  ::QuickAdvance( TemplateFieldTrack<Real_v>&       y_posvel,         // INOUT
                  const Real_v  dydx[],
                        Real_v  hstep,       // In
                  // Real_v& dchord_step,
                        Real_v& dyerr_pos_sq,
                        Real_v& dyerr_mom_rel_sq )
{
  // typedef typename Real_v Real_v;
//  typedef typename Mask<Real_v>      Bool_v;
  // Real_v dyerr_pos_sq, dyerr_mom_rel_sq;
  Real_v yerr_vec[TemplateFieldTrack<Real_v>::ncompSVEC],
           yarrin  [TemplateFieldTrack<Real_v>::ncompSVEC],
           yarrout [TemplateFieldTrack<Real_v>::ncompSVEC];
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
  y_posvel.LoadFromArray( yarrout, Nvar ); // fNoIntegrationVariables );
  y_posvel.SetCurveLength( s_start + hstep );

#ifdef  GUDEBUG_FIELD
  if(fVerboseLevel>2)
  {
    cout << "G4MagIntDrv: Quick Advance" << std::endl;
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
  dyerr_len = Sqrt( dyerr_len_sq );   // vecgeom::
  dyerr_len_sq /= epsilon ;

  // Set suggested new step
  hstep= ComputeNewStepSize( dyerr_len, hstep);
#endif

  return true;
}
 *********************************/

// --------------------------------------------------------------------------
#ifdef QUICK_ADV_ARRAY_IN_AND_OUT
template <class Real_v, class T_Stepper, unsigned int Nvar>
typename Mask<Real_v> RollingIntegrationDriver<Real_v, T_Stepper, Nvar>::QuickAdvance(Real_v yarrin[], // In
                                                                                     const Real_v dydx[],
                                                                                     Real_v hstep, // In
                                                                                     Real_v yarrout[],
                                                                                     Real_v &dchord_step,
                                                                                     Real_v &dyerr) // In length
{
  std::cerr << "ERROR in RollingIntegrationDriver::QuickAdvance()" << std::endl;
  std::cerr << "      Method is not yet implemented." << std::endl;

  //            FatalException, "Not yet implemented.");
  dyerr = dchord_step = hstep * yarrin[0] * dydx[0];
  yarrout[0]          = yarrin[0];
  exit(1);
}
#endif

/*****************************************
// --------------------------------------------------------------------------

//  This method computes new step sizes - but does not limit changes to
//  within  certain factors
//
template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
Real_v RollingIntegrationDriver< T_Stepper, Nvar>::
   ComputeNewStepSize(
      Real_v errMaxNorm,   // max error  (normalised)
      Real_v hStepCurrent) // current step size
{
  using Bool_v = vecCore::Mask_v<Real_v>;
  const int powerShrink = GetPowerShrink();
  const int powerGrow = GetPowerGrow();

  Bool_v goodStep = (errMaxNorm <= 1.0);
  Real_v powerUse = vecCore::Blend(goodStep, powerShrink, powerGrow);
  Real_v stretch  = GetSafetyFactor() * vecgeom::Pow(errMaxNorm, powerUse);
  Real_v hNew     = stretch * hStepCurrent;
  return hNew;
}

// ---------------------------------------------------------------------------

// This method computes new step sizes limiting changes within certain factors
//
// It shares its logic with AccurateAdvance.
// They are kept separate currently for optimisation.
//
template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
Real_v RollingIntegrationDriver<T_Stepper, Nvar>::
   ComputeNewStepSize_WithinLimits(
      Real_v errMaxNorm,   // max error  (normalised)
      Real_v hStepCurrent) // current step size
{
  using Bool_v = vecCore::Mask_v<Real_v>;

  Bool_v goodStep = (errMaxNorm <= 1.0);
  Real_v powerUse = vecCore::Blend(goodStep, fPowerShrink, fPowerGrow);

  Real_v stretch = fSafetyFactor * vecgeom::Pow(errMaxNorm, powerUse);

  Real_v stemp;
  stemp   = vecCore::math::Max(stretch, fMaxSteppingDecrease);
  stretch = vecCore::math::Min(stemp, fMaxSteppingIncrease);

  Real_v hNew = stretch * hStepCurrent;

#if 0
    // Compute size of next Step for a failed step
    if (errMaxNorm > 1.0 )
    {
      // Step failed; compute the size of retrial Step.
      hnew = fSafetyFactor * hstepCurrent * Math::Pow(errMaxNorm,fPowerShrink) ;

      hnew = std::min( hnew, fMaxSteppingDecrease * hstepCurrent );
                           // reduce stepsize, but no more
                           // than this factor (value= 1/10)
      }
    }
    else
    {
      // Compute size of next Step for a successful step
      if (errMaxNorm > fErrcon)
       { hnew = fSafetyFactor * hstepCurrent * Math::Pow(errMaxNorm,fPowerGrow); }
      else  // No more than a factor of 5 increase
       { hnew = fMaxSteppingIncrease * hstepCurrent; }
    }
#endif

  return hNew;
}

// ---------------------------------------------------------------------------
template <class T_Stepper, unsigned int Nvar>
void RollingIntegrationDriver<T_Stepper, Nvar>::ReportInvalidStepInputs(double hStepArr[], int nTracks)
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
******************************************/

// ------------------------------------------------------------

// ####################  Testing method(s) ####################

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
bool RollingIntegrationDriver<T_Stepper, Nvar>::TestInitializeLanes() // int numTracks)
{
  using std::cout;
  using std::cerr;
  using std::endl;
  using vecCore::Get;

  bool allOk = true;

  cout << " TestInitializeLanes() called. " << endl;

  constexpr int numTracks = 8;
  FieldTrack yInput[numTracks];
  // bool       badStepSize[numTracks];
  // double xStart[numTracks] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
  double chargeArr[numTracks] = {1.0, -2.0, 3.0, -4.0, 5.1, -6.1, 7.1, -8.1};
  // ------------------------ Test case 1 ------------------------
  Real_v yWorkLanes[Nvar], hStepLane, chargeLane, sBegin;
  double hStep[numTracks] = {0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1};
  //     ======
  int numBadSize = -1;

  cout << "    calling CreateInput() " << endl;
  CreateInput(yInput, numTracks);

  cout << "    test:  chargeArr [0] = " << chargeArr[0] << " ,  [1] = " << chargeArr[1] << ",  [2] = " << chargeArr[2]
       << ",  [3] = " << chargeArr[3] << endl;
  cout << "    calling InitializeLanes() " << endl;
  int nFilled1 = -1;
  int numNext1 = InitializeLanes(yInput, hStep, chargeArr, /* xStart, */ numTracks, // badStepSize,
                                 yWorkLanes, hStepLane, chargeLane, sBegin, nFilled1, numBadSize);
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
  int nextLocation2 = InitializeLanes(yInput, hStep2, chargeArr, /* xStart, */ numTracks, // badStepSize,
                                      yWorkLanes2, hStepLane, chargeLane, sBegin, nFilled2, numBadSize);

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
  ok1 = CheckOutput(yWorkLanes2, lane = 0, slot = 2, testName, varNameOutput);
  ok2 = CheckOutput(yWorkLanes2, lane = 1, slot = 4, testName, varNameOutput);
  ok3 = CheckOutput(yWorkLanes2, lane = 2, slot = 6, testName, varNameOutput);
  allOk = allOk && ok1 && ok2 && ok3;

  return allOk;
}

template <class T_Stepper, unsigned int Nvar>
void RollingIntegrationDriver<T_Stepper, Nvar>::AccurateAdvance(const FieldTrack yInput[], const double hstep[],
                                                               const double charge[],
                                                               double epsilon, // Can be scalar or varying
                                                               FieldTrack yOutput[], int nTracks, bool stillOK[]) const
{
  AccurateAdvance<geant::Double_v>(yInput, hstep, charge,
                                   epsilon, // Can be scalar or varying
                                   yOutput, stillOK, nTracks);
}

#ifdef EXTEND_SINGLE
template <class T_Stepper, unsigned int Nvar>
void RollingIntegrationDriver<T_Stepper, Nvar>::AccurateAdvance(const FieldTrack &yInput, const double hstep,
                                                               const double charge, double epsilon, FieldTrack &yOutput,
                                                               bool succeeded) const
{
  AccurateAdvance<double>(&yInput, &hstep, &charge,
                          epsilon, // Can be scalar or varying
                          &yOutput, &succeeded, 1);
}
#endif

#endif /* RollingIntegrationDriver_Def */
