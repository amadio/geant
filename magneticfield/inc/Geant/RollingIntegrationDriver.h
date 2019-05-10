//
// Driver for explicit Runge-Kutta methods
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and ensures
//  that the errors are within acceptable bounds.
// When multiple tracks are integrated, allows all tracks to continue
//   until all have completed a successful step and
//   at least one has completed its full integration interval (not just
//   the current step).
//
// History:
//
// Second Integration driver -  J. Apostolakis Nov / Dec 2018
//     Uses / extends idea of 'KeepStepping', first implemented by
//     Ananya Feb/March 2016
//
// Includes changes from the version of SimpleIntegrationDriver adapted for modularisation
//     April 2019
//
// First templated version:  Ananya, Feb/March 2016
//     ( commit 95e1316bcc156a04c876d6ea0fc9e60a15eeac4f )
//
// - Contributors: J.Apostolakis                    Nov 2018 - Jan 2019
// --------------------------------------------------------------------

#ifndef RollingIntegrationDriver_Def
#define RollingIntegrationDriver_Def

// #include "Geant/TemplateFieldTrack.h"
#include "base/AlignedBase.h"
#include "base/Vector.h"

#include "Geant/math_wrappers.h"

#include "Geant/FieldTrack.h"

// #include "TemplateVScalarIntegrationStepper.h"
// #include "IntegrationStepper.h"


// Adding because adding scalar stepper for new constructor (KeepStepping)
// #include "Geant/VScalarIntegrationStepper.h"

// Adding to send in scalar driver to deal with 1/2 remaining lanes
// #include "Geant/ScalarIntegrationDriver.h"
// #include "Geant/ScalarFieldTrack.h"

#include "Geant/FlexIntegrationDriver.h"

#include "Geant/BaseRkIntegrationDriver.h"
#include "Geant/AuxVecMethods.h"

// ----  Flags for debugging and diagnostics only - off in benchmarks!
#define  DRIVER_PRINT_PROGRESS  1
#define  DRIVER_DIAGNOSTICS     1

#include "ErrorEstimatorSixVec.h"
    
// #include "Geant/FormattedReporter.h"  // Direct include is not needed.
#include "Geant/PrintDriverProgress.h"

#include "Geant/VectorTypes.h" //  Defines geant::Double_v
#include "Geant/math_wrappers.h"

#define CHECK_ONE_LANE 1
//  Define to check a single lane

#define CONST_DEBUG    1
//  Define to turn 'partDebug' into compile time constant

#ifdef CHECK_ONE_LANE
#include "IntegrationDriverConstants.h"
#endif

// --------------------------------------------------------------
template <class T_Stepper, unsigned int Nvar>
class RollingIntegrationDriver :
   // public FlexIntegrationDriver,
   public BaseRkIntegrationDriver<T_Stepper, Nvar> // ,
   // public vecgeom::AlignedBase
{
public:
  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

  RollingIntegrationDriver(double     hminimum, // same
                           T_Stepper *pStepper,
                           double     epsRelativeMax, // = 1.0e-5,
                           int        numberOfComponents = 6 );

  virtual ~RollingIntegrationDriver();

  virtual void AccurateAdvance(const FieldTrack yInput[], const double hstep[], const double charge[], // double epsilon, => moved to c-tor
                               FieldTrack yOutput[], int nTracks, bool succeeded[]) const override final;

#ifdef EXTEND_SINGLE
  virtual void AccurateAdvance(const FieldTrack &yInput, const double hstep, const double charge, // double epsilon,
                               FieldTrack &yOutput, bool succeeded) const override final;
#endif

  // Implemented in terms of the following templated function:
  template <class Real_v>
  void AccurateAdvance(const FieldTrack yInput[], const double hstep[], const double charge[], // double epsilon,
                       FieldTrack yOutput[], bool succeeded[], int nTracks) const;
  // Drive Runge-Kutta integration of ODE for several tracks (ntracks)
  // with starting values yInput, from current 's'=0 to s=h with variable
  // stepsize to control error, so that it is bounded by the relative
  // accuracy eps.  On output yOutput is value at end of interval.
  // The concept is similar to the odeint routine from NRC 2nd edition p.721

  // Auxiliary methods
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
  void KeepStepping(const Real_v   yStart[], //  [N]
                    Real_v         dydx[],   //  [N]
                    const Real_v   charge,
                    bool           lastBatch,  // There is no more work after this -- so keep going at end (or use Sequential)
                    Real_v       & xVal, // InOut
                    Real_v         xEnd, // Input only - End of Interval                   
                    Real_v         htry,
                    // double         epsilon, // req. relative accuracy - Was const Real_v
                    Real_v         yEnd[],  // [N]
                    Real_v       & hdid,
              //    Real_v       & hLast,   // Out - last attempted step
                    Real_v       & hnext) const;
  // Integrate over the length htry 'together', while each lane succeeds,
  //   until the ends of at least one interval (xEnd) or failure (underflow). 
  // The initial length 'htry', can be reduced or increased for
  //   subsequent integration steps.
  // 
  // In this version it will return as soon as one lane ends or fails.

  template <class Real_v>
  int InitializeLanes(const FieldTrack yInput[], const double hstep[], const double charge[], int nTracks,
                      int indexArr[], // [vecCore::VectorSize<Real_v>()] - Output
                      Real_v y[],     // [Nvar]        - Output
                      Real_v &hStepLane, Real_v &chargeLane, Real_v &startCurveLength,
                      int &numFilled,           // [Out]: number of filled lanes
                      bool     badStepSize[],   // Output - needed to reset values                      
                      int &numBadSize) const;
  // Load 'work' into array y, indices of lanes into indexArr, etc
  // Return array index of 'next' track - i.e. first track not yet loaded

  template <class Real_v>
  bool InsertNewTrack(const FieldTrack yInput[], const double hstep[], const double charge[], const int slot,
                      int    & trackNextInput, bool succeeded[], Real_v y[], Real_v &hStepLane, Real_v &chargeLane,
                      Real_v & startCurveLength,
                      int      indexArr[],
                      int      nTracks,
                      bool     badStepSize[],  
                      int    & numBadSize
     ) const;

  template <class Real_v>
     void StoreOutput(const Real_v yEnd[],     // End state  ( indOut )
                      const Real_v x,          // End curve length
                      int          currIndex,  // Index of lane (within Real_v )  - make it unsigned !
                      FieldTrack   yOutput[],
                      int          indOut,     // location in Output arrays ( 0 <= value < nTracks )
                      const double hstep[],    // Don't store if value hstep[ indOut ] <= 0.0
                      bool         succeeded[],
                      int          nTracks) const;
  // Called whenever a lane is finished, to store output of integration
  // Input arguments as above.  
  //    - Note: will not store if hstep[indOut] <= 0.0
  // Stores values into: 
  //    - yOutput:  end position, momentum and curve length
  //    - succeeded[nTracks]  success flag

  void CheckParameters()     { BaseRkIntegrationDriver<T_Stepper, Nvar>::CheckParameters(); }
  int  IncrementStepperCalls() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::IncrementStepperCalls() ; }
  int  IncrementNumberSteps()  { return ++fNoTotalSteps; }

  template <class Real_v>  
    Real_v ComputeNewStepLengthWithinLimits2( Real_v errMaxSq, Real_v hStepCurrent) const
      { return BaseRkIntegrationDriver<T_Stepper, Nvar>::
                       ComputeNewStepLengthWithinLimits2( errMaxSq, hStepCurrent ); }

  template <class Real_v>  
    Real_v ComputeNewStepLengthWithinLimits3( Real_v    errMaxSq,
                                              Real_v    hStepCurrent,
                                              vecCore::Mask_v<Real_v> isNeeded,  // Bool_v
                                              Real_v  & stretchFactor // OUT: stretch (if needed only)
       ) const
      { return BaseRkIntegrationDriver<T_Stepper, Nvar>::
            ComputeNewStepLengthWithinLimits3( errMaxSq, hStepCurrent, isNeeded, stretchFactor ); }
     
   // Calculated new step limit (and stretch factor) - but only if 'isNeeded' is true && (hStepCurrent != 0.)xo
   
  // Setting parameters ( few / none now )

  // Compute dependent parameters

  // using Bool_v = vecCore::Mask_v<Real_v>;
#ifdef ERRCON_INBASE   
  void ComputeAndSetErrcon() { BaseRkIntegrationDriver<T_Stepper, Nvar>::ComputeAndSetErrcon(); }
  // double GetErrcon()      const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetErrcon(); }
  double GetErrcon()   const {
     double val= BaseRkIntegrationDriver<T_Stepper,Nvar>::GetErrcon();
     std::cout << " <--- BaseRkIntDriver gives errcon = " << fErrcon << " ---> "; 
     assert( val > 0.0 );
     return val;
  }  
#else
  void ComputeAndSetErrcon() {
      fErrcon       = Math::Pow(fMaxSteppingIncrease / fSafetyFactor, 1.0 / GetPowerGrow() );
      fShrinkThresh = Math::Pow(fMaxSteppingDecrease / fSafetyFactor, 1.0 / GetPowerShrink() );
      
      std::cout << "SimpleIntegrationDriverComputAndSetErrcon():  fErrcon = " << fErrcon
                << "  from:  maxStepIncrease =  " << fMaxSteppingIncrease
                << "  fSafetyFactor = " << fSafetyFactor
                << "  power-grow =  " << GetPowerGrow() << std::endl;
      
      assert( fErrcon > 0.0 ); 
  }
  double GetErrcon()       const { assert( fErrcon > 0.0 ); return  fErrcon; }
  double GetShrinkThresh() const { assert( fShrinkThresh > 0.0); return fShrinkThresh; }
   
private: 
  double fErrcon= 0.0;
  double fShrinkThresh= 0.0; 
#endif   

   
protected:
  template <class Real_v>
  void Report1(const vecCore::Mask_v<Real_v>   finished,   // Bool_v
               const vecCore::Mask_v<Real_v>   active,     // Bool_v
               const Real_v charge,
               const Real_v yEndStep[],
               const Real_v yErr[],
               const Real_v h  ) const
      { PrintDriverProgress::Report1<Real_v,Nvar>( finished, active, charge, yEndStep, yErr, h ); }
  
public:
  // Access parameters
  double GetPowerShrink()  const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fPowerShrink; }
  double GetPowerGrow()    const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::fPowerGrow; }
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

   // double GetMaxRelativeEpsilon() const { return fEpsilonRelMax; }

  inline unsigned int GetMaxNoSteps() const { return fMaxNoSteps; } 
  // inline unsigned int GetMaxNoSteps() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetMaxNoSteps() ; } // fMaxNoSteps; }  
  
  // Modifier methods
  void SetMaxNoSteps(int val) { BaseRkIntegrationDriver<T_Stepper,Nvar>::SetMaxNoSteps(val) ; }
  
  int IncrementStepperCalls()  { return BaseRkIntegrationDriver<T_Stepper,Nvar>::IncrementStepperCalls() ; }

  // const T_Stepper *GetStepper() const { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStepper(); }
  //       T_Stepper *GetStepper()       { return BaseRkIntegrationDriver<T_Stepper,Nvar>::GetStepper(); }
  
  const T_Stepper *GetStepper() const { return fpStepper; }   // Tried suppressing to check its use  2019.04.10 @ 17:20
        T_Stepper *GetStepper() { return fpStepper; }

  void CreateInput(FieldTrack yInput[], int nTracks)
  { BaseRkIntegrationDriver<T_Stepper,Nvar>::CreateInput( yInput, nTracks);  }
        
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

private:
  // Private copy constructor and assignment operator.

  RollingIntegrationDriver(const RollingIntegrationDriver &);
  // Copy constructor used to create Clone method

  RollingIntegrationDriver &operator=(const RollingIntegrationDriver &) = delete;

  template <class Real_v>
  void StoreAllValues( const Real_v                    yWork[], 
                       const Real_v                  & hValue,
                       const Real_v                  & errMaxSqValue,
                       Real_v                          yFinal[],
                       Real_v                        & hFinalStore,
                       Real_v                        & errMaxSqStore
     ) const;
   
  template <class Real_v>
  void StoreGoodValues(const Real_v                    yWork[],
                       const Real_v                  & hValue,
                       const Real_v                  & epsSqVal,
                       const vecCore::Mask_v<Real_v> & storeFlag,
                       Real_v                          yFinal[],
                       Real_v                        & hFinalStore,
                       Real_v                        & epsSqStore
     ) const;
  // Auxiliary method, to store results of selected ('good') lanes

private:
  // ---------------------------------------------------------------
  // DEPENDENT Objects
  T_Stepper                   *fpStepper;
  const  ErrorEstimatorSixVec  fErrorEstimator;

  // ---------------------------------------------------------------
  //  INVARIANTS

  // Moved to BaseRkIntegrationDriver
  //
  // double fMinimumStep; // same
  //    Minimum Step allowed in a Step (in absolute units)
  // static constexpr double fSmallestFraction = 1.0e-7; // Expected value: larger than 1e-12 to 5e-15;
  //    Smallest fraction of (existing) curve length - in relative units
  //    below this fraction the current step will be the last

  // const double fEpsilonRelMax; //  Maximum Relative Error in integration ==> Moved to FlexIntegrationDriver

  const double fHalfPowerShrink;
  unsigned long fMaxNoSteps= 100;  // Default - expect it to be overriden in constructor

  // ---------------------------------------------------------------
  // Compilation constants
#ifdef CONST_DEBUG  
  const bool partDebug  = true;                 // Enforce debugging output
#endif
  // const int ncompSVEC   = FieldTrack::NumCompFT; // expect 6, later 8, eventually up to 12
  const bool useOneStep = true;                  //  Algorithm selection - false for KeepStepping

  // ---------------------------------------------------------------
  //  INVARIANTS - continued 

  // const int  fNoIntegrationVariables;  // Number of Variables in integration
  // const int fMinNoVars; // Minimum number for TemplateFieldTrack<Real_v>
  // const int fNoVars;    // Full number of variable

  static constexpr int fMaxStepBase = 250;

  // static constexpr double fSafetyFactor= 0.9; // -> Fails to compile on clang 9.1 2017.12.05
  const double fSafetyFactor = 0.9; //     OK ...
  // const double fPowerShrink;        //  exponent for shrinking
  // const double fPowerGrow;          //  exponent for growth
  // Parameters used to grow and shrink trial stepsize.

  // double fSurfaceTolerance = 1.e-6;

  //  Stepsize can increase by no more than 5.0
  //           and decrease by no more than x10. = 0.1
  static constexpr double fMaxSteppingIncrease = 5.0;
  static constexpr double fMaxSteppingDecrease = 0.1;
  // Maximum stepsize increase/decrease factors.

  // static constexpr double fSmallestFraction = 1.0e-7; // Expected value: larger than 1e-12 to 5e-15;
  const double fSmallestFraction = 1.0e-7; // Expected value: larger than 1e-12 to 5e-15;
   
  // int fStatisticsVerboseLevel         = 0;
  mutable std::atomic<unsigned long> fStepperCalls;
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

  mutable unsigned long fNoTotalSteps = 0; // , fNoBadSteps = 0, fNoSmallSteps = 0, fNoInitialSmallSteps = 0;

  int fVerboseLevel; // Verbosity level for printing (debug, ..)
                     // Could be varied during tracking - to help identify issues

}; // End of class definition -- RollingIntegrationDriver

#ifdef DRIVER_PRINT_PROGRESS  
#include "Geant/FormattedReporter.h"
#endif

/***************************

template <class T_Stepper, unsigned int Nvar>
inline void RollingIntegrationDriver<T_Stepper, Nvar>::CheckParameters()
{
  constexpr double perMillion = 1.0e-6;
  using std::cerr;
  using std::endl;

  double checkPowerShrink = -1.0 / fpStepper->GetIntegratorOrder();

  double diffShrink = fPowerShrink - checkPowerShrink;
  if (std::fabs(diffShrink) // checkPowerShrink - fPowerShrink)
      >= perMillion * std::fabs(fPowerShrink)) {
    cerr << "RollingIntegrationDriver: ERROR in fPowerShrink" << std::endl;
    cerr << "    calculated = " << checkPowerShrink << "    pre-computed = " << fPowerShrink << "  diff= " << diffShrink
         << "  tolerance = " << perMillion * std::fabs(fPowerShrink) << endl;
    cerr << "  Order of integrator = " << fpStepper->GetIntegratorOrder() << endl;
    exit(1);
  }
  assert(std::fabs(checkPowerShrink - fPowerShrink) < perMillion * std::fabs(fPowerShrink));

  double checkPowerGrow = -1.0 / (1.0 + fpStepper->GetIntegratorOrder());
  assert(std::fabs(checkPowerGrow - fPowerGrow) < perMillion * std::fabs(fPowerGrow));

  if (std::fabs(checkPowerGrow - fPowerGrow) >= perMillion * std::fabs(fPowerGrow)) {
    std::cerr << "RollingIntegrationDriver: ERROR in fPowerGrow" << std::endl;
    exit(1);
  }

  if (fVerboseLevel)
    std::cout << "RollingIntegrationDriver::CheckParameters > Powers used: " << std::endl
              << "  shrink = " << fPowerShrink << "  grow = " << fPowerGrow << std::endl;
}
***********************/


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
RollingIntegrationDriver<T_Stepper, Nvar>::
      RollingIntegrationDriver(double      hMinimumStep,
                               T_Stepper * pStepper,
                               double      epsRelMax,
                               int         numComponents )
    :
      BaseRkIntegrationDriver<T_Stepper, Nvar>(hMinimumStep,
                                               pStepper,     // or pass pStepper->GetIntegratorOrder()  ??
                                               epsRelMax,
                                               numComponents,
                                               false),  // statistics Verbosity
      // fEpsilonRelMax( epsRelMax ),
      fErrorEstimator( epsRelMax, hMinimumStep ),
      fHalfPowerShrink( 0.5 * BaseRkIntegrationDriver<T_Stepper, Nvar>::GetPowerShrink() )
{
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
  assert(pStepper != nullptr);
  assert(Nvar <= (unsigned int)numComponents /* was ncompSVEC */ ); // Ensure that arrays are large enough for Integr.

  fpStepper = pStepper;

  // ComputeAndSetErrcon();
  
  // fMaxNoSteps = fMaxStepBase / fpStepper->GetIntegratorOrder();
  SetMaxNoSteps( fMaxStepBase / fpStepper->GetIntegratorOrder() );

  ComputeAndSetErrcon();
  assert( GetErrcon() > 0.0 );
  
  CheckParameters();

#ifdef GUDEBUG_FIELD
  fVerboseLevel = 2;
#endif

  if (fVerboseLevel) {
    std::cout << "RollingIntegrationDriver/RiD:ctor> Stepper Order= " << pStepper->GetIntegratorOrder() << " > Powers used: "
              << " shrink = " << GetPowerShrink() << "  grow = " << GetPowerGrow()
              << " and with  max-relative-error = " << epsRelMax       
              << std::endl;
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
    const RollingIntegrationDriver<T_Stepper, Nvar> &right)
    :
   BaseRkIntegrationDriver<T_Stepper, Nvar>::BaseRkIntegrationDriver( right ) 
   // , fStatisticsVerboseLevel(right.fStatisticsVerboseLevel)
{
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
   
  // const T_Stepper *protStepper = right.fpStepper;
  // fpStepper                    = protStepper->Clone();
  // ==> This should happen in the base class .....
   
  ComputeAndSetErrcon();
  assert( GetErrcon() > 0.0 );
  
  SetMaxNoSteps( std::max( fMaxStepBase / GetStepperOrder(), 10 ) );  // fpStepper->GetIntegratorOrder() );

  if(GetVerboseLevel() > 0) {
    std::cout << "RollingIntegrationDriver copy Constructor. " << std::endl;
  }
}

// ---------------------------------------------------------

//  Destructor
template <class T_Stepper, unsigned int Nvar>
RollingIntegrationDriver<T_Stepper, Nvar>::
   ~RollingIntegrationDriver()
   //======================
{
   if (GetStatisticsVerboseLevel() > 1) {
      // PrintStatisticsReport();
   }
}

// ---------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void RollingIntegrationDriver<T_Stepper, Nvar>::
   KeepStepping(const Real_v   yStart[], //  [N]
                      Real_v   dydx[],
                const Real_v   charge,
                bool           lastBatch,    // Last batch of tracks - finish it before returning !
                Real_v       & x,    // InOut
                Real_v         xEnd, // Input only - End of Interval
                Real_v         htry,
                // double        eps_rel_max,  // ==> Now a class constant (embedded in Error Estimator)
                // const Real_v  eps_rel_max,  //     Older choice - different value per lane
                Real_v         yFinal[], // Out-values
                Real_v       & hdid,    // Out - total achieved length
//              Real_v       & hlast,   // Out - length of last (attempted) step
                Real_v       & hnext    // Out - proposed next integration length
   )
   const
// This version:  J. Apostolakis,  13 November 2017.
//   Lanes are integrated until
//     -  finish their integration ( reach xEnd for the lane ), or
//     -  fail due to step-size underlow.
//   or else
//     -  a maximum number of trial iterations is reached,
//  (No minimum step size is used.)
//
// Refined method, suitable for cases where lanes have different
//   lenghts (in angle), or difficulty to succeed.
// -------------------------------------------------------------------------
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
  using FormattedReporter::ReportRowOfDoublesIf;
  using FormattedReporter::ReportOneLane;
  using PrintDriverProgress::ReportConditionLanes;

  if (partDebug) { cout << "\n" << endl; }

  Real_v yerr[Nvar];
  Real_v yStepStart[Nvar];
  Real_v yStepEnd[Nvar];    // Could reuse yFinal for this ??  2018.01.08
  Real_v hStep = htry; // Set stepsize to the initial trial value
  
  const Real_v xStart= x;    // Remember

  // Store values for last good step - for use when current step is not good !
  Real_v hStepLastGood(0.0);
  Real_v errmaxSqLastGood(-1.0);
  
  static int tot_no_trials = 0; // Should be thread_local - or suppressed. Just statistics
  const int max_trials     = 100;

  Bool_v finishedLane   = (htry <= 0.); //  Allows h <=0 as signal lane is empty.
  Bool_v progressedLane(false);
  
  Real_v hFinal = htry, errmax_sqFinal = Real_v(0.); // xFinal, hdidFinal,
  Real_v hnew;
  Real_v errmaxSqFallThru(0.0);  

  Bool_v goodStep(false);
  Bool_v UnderflowOccurred(false);
  
  // Bool_v keepGoing(true);
  bool   alreadySomeDone= false;
  //  bool   doneExtraIteration= false;
  
  // for (iter = 0; iter < max_trials; iter++)
  unsigned int iter = 0;
  int itersLeft     = max_trials;
  // ReportManyRowsOfDoubles( "yStart",  yStart, Nvar );
  Real_v magmomInit_sq = yStart[3] * yStart[3] + yStart[4] * yStart[4] + yStart[5] * yStart[5];

#ifdef CHECK_ONE_LANE
  const int trackToPrint = IntegrationDriverConstants::GetInstance()->GetTrackToCheck();
#ifndef CONST_DEBUG
  partDebug = ( laneToCheck != -1 );
#endif
  // bool  printNow = (laneToCheck >= 0) && vecCore::Get( htry, laneToCheck ) > 0.0 ;       
  if( laneToCheck != -1 ) {
     cout << "* RID:AccAdv Global Vars> trackToPrint = " << trackToPrint << " l2c / laneToCheck = " << laneToCheck;
     cout << "  Args>  h[l2c] = " << vecCore::Get( htry, laneToCheck );     
     cout << endl;
  }
#endif
  if( partDebug ) { std::cout << "RID: KeepStepping called." << std::endl; }
       
  // ErrorEstimatorSixVec fErrorEstimator( eps_rel_max, GetMinimumStep() );

  for( unsigned int i=0; i<Nvar; i++)
    yStepStart[i] = yStart[i];

  // Bool_v lastStepWasGood(false);
  Bool_v integrationDone(false);
  
  do {
    Bool_v Active = !finishedLane;
    Real_v errmax_sq = 0.0;

#ifdef CHECK_ONE_LANE    
    bool  printLane = (laneToCheck >= 0) && vecCore::Get( Active, laneToCheck );
#endif        

    vecCore::MaskedAssign(hStep, finishedLane, Real_v(0.0)); // Set hStep = 0.0 for finishedLane lanes

    bool anyProgressLastStep = ! vecCore::MaskEmpty(goodStep) ;
    if( iter > 0 ) // && ! vecCore::MaskEmpty(goodStep) )  // If none moved, don't calculate derivatives
    { 
       fpStepper ->RightHandSideInl( yStepStart, charge, dydx ); // TODO: change to inline
       //---------****************-----------------------
       cout << "RID> Recalculated   dy/ds ";
       if( !anyProgressLastStep ) { cout << " (likely unnecessary) "; }
       else                       { cout << " (believed needed)    "; }
    } else {
       if( iter > 0 ) { cout << "RID> Kept values of dy/ds from previous step "; }
    }
    cout << " - iter = " << iter << " anyProgressLastStep = " << anyProgressLastStep << std::endl;       
    ReportManyRowsOfDoubles("Current X/P",  yStepStart, 6 );
    ReportRowOfDoubles("Charge",  charge, 6 );       
    ReportManyRowsOfDoubles("d/ds [X,P] ", dydx, 6 );

    itersLeft--;
    iter++;
    if (partDebug) cout << " KeepStepping - iteration = " << iter << "             "
                        << " ( total iterations = " <<  (tot_no_trials + iter) << " ) "
                        << endl;
    
    fpStepper->StepWithErrorEstimate(yStepStart, dydx, charge, hStep, yStepEnd, yerr); // CAREFUL -> changes for others ?
    //******************************
    
    IncrementStepperCalls();    //  fStepperCalls++;

#ifdef DRIVER_PRINT_PROGRESS
    bool DebugEachIteration = false;
    if (partDebug && DebugEachIteration) { Report1( finishedLane, Active, charge, yStepEnd, yerr, hStep ); }
#endif

#ifdef CHECK_ONE_LANE
    Real_v errpos_sq = 0.0; // square of displacement error
    Real_v errmom_sq = 0.0; // square of momentum vector difference
    Real_v epsPosition = 0.0; 
    errmax_sq = fErrorEstimator.EstimateError( yerr, hStep, magmomInit_sq, epsPosition, errpos_sq, errmom_sq );
#else
    errmax_sq = fErrorEstimator.EstimateError( yerr, hStep, magmomInit_sq ); // Older ver:  yStepEnd );
#endif

    goodStep = Active && (errmax_sq <= 1.0);
    progressedLane =  progressedLane || goodStep;

    vecCore::MaskedAssign( hStepLastGood,    goodStep, hStep );     // if( goodStep[i] ) hStepLastGood[i] = hStep[i]; 
    vecCore::MaskedAssign( errmaxSqLastGood, goodStep, errmax_sq ); // if( goodStep[i] ) errmaxSqLastGood[i] = errmax_sq[i]; 
    
    // Update hdid for successful steps
    Real_v hAdvance= vecCore::Blend( goodStep, hStep, Real_v(0.0) );
    hdid += hAdvance;

    // Real_v xNext = x + hStep; //   Updates both good and bad steps -- else it should be  = x + hAdvance;
    Real_v xNext = x + hAdvance; //   Updates only good steps --  2019.01.07

    // fSmallestFraction    
    Bool_v endingIntegration  = ( xEnd - xNext < geant::units::perMillion * htry ) ;  // Care for very small diff
    Bool_v laneDone = ( endingIntegration | finishedLane);      // WAS = (goodStep | finished);

    integrationDone =  integrationDone | endingIntegration;
  
    if( partDebug ) {
       cout << " ====================================================================" << endl;
       ReportRowOfDoubles("charge ***>", charge );           
       ReportRowOfDoubles("x:   Start  =", x );
//       ReportRowOfDoubles("xEnd: ok/bad=", x + hStep );
       ReportRowOfDoubles("xNext: End  =", xNext );
       ReportRowOfDoubles("hStep/tried =", hStep );
       ReportRowOfDoubles("hAdvance/got=", hAdvance );
       ReportRowOfDoubles("hdid (total)=", hdid );
       // ReportRowOfDoubles("errmax^2 =", errmax_sq );       // ADD soon
       ReportRowOfBools<Real_v>("goodStep", goodStep);
       ReportRowOfBools<Real_v>("progressedLane", progressedLane);       
       ReportRowOfBools<Real_v>("endingIntegration", endingIntegration);
       ReportRowOfBools<Real_v>("laneDone", laneDone);              
       cout << " ====================================================================" << endl;
    }
    
    bool  allDone =    vecCore::MaskFull(laneDone);
    bool  enoughFinished = allDone;  // Not final value, except to experiment
    bool  someDone =   allDone; // Ver A.1  Original:  ! vecCore::MaskEmpty(laneDone);
    if( !allDone ) {
       someDone = ! vecCore::MaskEmpty(laneDone); // Ver A.1
       // enoughFinished = someDone;             // First Try
       // int   countDone = countLanesTrue( laneDone );
       enoughFinished = allDone || (someDone && alreadySomeDone);
    }
    // bool enoughFinished = someDone; // First Try    
    // bool enoughFinished = allDone || (someDone && allreadySomeDone);
    //
    // Ideas for refinement of criterion to finish:
    //  v2   - take exactly one more step if some had finished             --> Current    
    //  v2.1 - take exactly one more step if only one finished last iteration

    bool mustFinishAllStill = lastBatch && ! allDone;  // For now cannot divert single/few lane into serial code
                                                  // So must do all the work here!

    std::cout << "Loop end condition: " << " enoughFinished= " << enoughFinished
              << " composed using : " 
              << " someDone= " <<  someDone
              << " allDone= " << allDone
              << " alreadySomeDone = " << alreadySomeDone << " (from last iteration) "
              << " mustFinishAllStill = " << mustFinishAllStill       
              << std::endl;
    
    // alreadySomeDone = someDone;    // ---> Moved to just before end of while loop
    
    finishedLane = laneDone;
    
#ifdef CHECK_ONE_LANE
    // Debugging one lane (at a time)  -------------  2019.04.17
    if( printLane ) { // if ( (laneToCheck >= 0) && vecCore::Get( Active, laneToCheck ) ) {
       // const int trackToPrint = IntegrationDriverConstants::GetInstance()->GetTrackToCheck();
       // bool    allDone =    vecCore::MaskFull(laneDone);
       ReportOneLane ( hStep, x, epsPosition, errpos_sq, errmom_sq, errmax_sq, laneDone,
                       allDone, iter, tot_no_trials, laneToCheck, trackToPrint,
                       "RollingID" );
       std::cout << " Track to check " << trackToPrint << " laneToCheck = " << laneToCheck << std::endl;
       ReportRowOfSquareRoots("ErrPos", errpos_sq );
       ReportRowOfSquareRoots("ErrMom", errmom_sq );
       ReportRowOfSquareRoots("ErrMax", errmax_sq );
       
       if( 1 ) {
          std::cout << "RID: Status after stepper call ----------------------------------------------" << std::endl;
          ReportManyRowsOfDoubles("RealStartX/P", yStart, 6 );          
          FormattedReporter::FullReport(yStepStart, charge, dydx, hStep, yStepEnd, yerr, errmax_sq, Active, goodStep);
          // ReportManyRowsOfDoubles("err-p/xyz", &yerr[3], 3 );
          // ReportRowOfSquareRoots("|err-p|", yerr[3]*yerr[3] + yerr[4]*yerr[4] + yerr[5]*yerr[5] );       
          // ReportRowOfDoubles("up = SumErr^2", sumerr_sq );
          // ReportRowOfDoubles("dwn= magMom^2+e", magmom_sq + tinyValue );
          // ReportRowOfDoubles("mul:1/e_vel^2", invEpsilonRelSq );
       }
    }
    // End debug code                 -------------  2019.04.17
#endif 
    
    Active   = !finishedLane;
    x = xNext;  // Move the starting point of the next step
    
#ifdef DRIVER_PRINT_PROGRESS
    if (partDebug) {
      ReportRowOfBools<Real_v>("goodStep", goodStep);
      ReportRowOfBools<Real_v>("laneDone", laneDone);
      ReportRowOfBools<Real_v>("finished(upd)", finishedLane);
      ReportRowOfBools<Real_v>("Active(upd)", Active);
      
      ReportRowOfDoubles("x: NOT Updated =", x );
    }
#endif

    if( mustFinishAllStill ){
       std::cout << "RollingID> Now in 'FinishAll' mode ******** " << std::endl;
    }
    
    if ( enoughFinished && ! mustFinishAllStill ) 
    {  // Idea 1.5
       if (partDebug) cout << "RID: Enough Finished (Store and Report OFF ) - just breaking." << endl;
       // StoreGoodValues(yStepEnd, hStep, errmax_sq, goodStep, yFinal, hFinal, errmax_sqFinal);
       //   Why store here too ??? - store only after loop exit !     2019.01.07   ( TBC )
       errmaxSqFallThru = errmax_sq;

       // std::cout << "** Exiting loop - Report of status -------------------------------------- " << std::endl;
       // FormattedReporter::FullReport(yStepStart, charge, dydx, hStep, yStepEnd, yerr, errmax_sq, Active, goodStep);

       break;
    }

    // Refill StepStart using values from StepEnd -- for successful lanes
    for( unsigned int i=0; i<Nvar; i++)
       vecCore::MaskedAssign( yStepStart[i], goodStep, yStepEnd[i] );

    Real_v hReduced = ComputeNewStepLengthWithinLimits2( errmax_sq, hStep );
    //  Expect all (or nearly all) have work -- so using vecCore::Pow on all lanes makes sense.  Else ...

    hnew = vecCore::Blend(finishedLane, Real_v(0.0), hReduced);    

    Real_v xnew        = x + hnew;

    Bool_v stepSizeUnderflow = Active && (xnew == x);

#ifdef CHECK_ONE_LANE    
    if( printLane ) {
      std::cout << "RID: After chance to break. (Not enoughFinished.) Remaining lanes represent work: " << std::endl;
      ReportRowOfBools<Real_v>("laneDone", laneDone);
      ReportRowOfBools<Real_v>("Active", Active);
      // ReportRowOfDoubles("powerShrink=", 2*fHalfPowerShrink );
      // std::cout << " safetyFactor = " << fSafetyFactor << std::endl;
      ReportRowOfDoubles("errMaxSq=", errmax_sq );
      ReportRowOfSquareRoots("errMax=", errmax_sq );
      // ReportRowOfDoubles("errPower=", errPower );  // Intermediate result. Subsumed in ComputeNewStepLengthWithinLimits2
      ReportRowOfDoublesIf<Real_v>("hReduced/ifActive", hReduced, Active);      
      // ReportRowOfDoubles("hReduced=", hReduced );
      ReportRowOfDoubles("hnew=", hnew );
      ReportRowOfDoubles("xnew=", xnew );
      if( ! vecCore::MaskEmpty(stepSizeUnderflow) ) {     
         ReportRowOfBools<Real_v>("underflow", stepSizeUnderflow);
         std::cout << "** WARNING: Underflow Detected in RollingIntegrationDrive !!" << std::endl;
         std::cerr << "RollingIntegrationDrivero> Underflow Detected !!" << std::endl;         
      }
    }
#endif
    
    finishedLane = finishedLane || stepSizeUnderflow;
    UnderflowOccurred = UnderflowOccurred || stepSizeUnderflow;

#ifdef DRIVER_DIAGNOSTICS  
    if (!vecCore::MaskEmpty(stepSizeUnderflow)) {
       Real_v xnewB       = x + hnew;     
       int numUnder= countMaskTrue<Real_v>( stepSizeUnderflow );
       std::cerr << "WARNING > In iteration loop: found " << numUnder
                 << "underflow lanes in this step." << std::endl;
       ReportConditionLanes(stepSizeUnderflow, x, xnewB, hStep, htry);
  }
#endif
    
    hStep = Min( hnew, xEnd - x ); // Ensure not to go past end

    // Interim Report of problem lanes -- For debugging only !!
    Bool_v problemLanes = stepSizeUnderflow;
    if (!vecCore::MaskEmpty(problemLanes)) {
      std::cerr << "GVIntegratorDriver::OneStep:" << std::endl
                << "  Stepsize underflow in Stepper ( Report 1 - in loop )" << std::endl;
      Bool_v problemLanes = stepSizeUnderflow && !finishedLane;

      ReportConditionLanes(problemLanes, x, xnew, hStep, htry);
    }
    
    // Refinement of criterion v2 - take one more step if only one is finished ... ?
    // if( someDone && ! alreadySomeDone ) {  alreadySomeDone = true; }
    alreadySomeDone = someDone;

  } while (itersLeft > 0 && (!vecCore::MaskFull(finishedLane))
           );
  // Only ensure that
  //   - iteration don't exceed max
  //   - some work remains to be done
  // An early exit from the 'break' condition is the 'normal' way out of the loop.
  
  tot_no_trials += iter;

  // Store exactly one time, here!

  // Original code - from Simple
  // StoreGoodValues(yStepEnd, hStep, errmaxSqFallThru, finishedLane, yFinal, hFinal, errmax_sqFinal);
  // Q> - Can shortcut this, but using yFinal[] instead of yStepEnd[] in code above.  Worthwhile or bad ? TODO

  // Try 2:   2019.05.07
  //   Store unconditionally - intermediate state must be starting point for next trial step !!
  // StoreGoodValues(yStepEnd, hStep, errmaxSqFallThru, Bool_v(true), yFinal, hFinal, errmax_sqFinal);
  //   WRONG: it also stores the result of failed steps.
  
  // Try 3:   2019.05.08
  // cout << " Store Good Values calls - version 3:   2019.05.08" << std::endl;  
  // StoreGoodValues(yStepEnd, hStep, errmaxSqFallThru, progressedLane, yFinal, hFinal, errmax_sqFinal);

#if 1
  // Try 4:   2019.05.09  11:20
  cout << " Store Good Values calls - version 4:   2019.05.09  11:20" << std::endl;
  StoreGoodValues(yStepEnd,   hStep, errmaxSqFallThru, progressedLane &&  goodStep, yFinal, hFinal, errmax_sqFinal);
  StoreGoodValues(yStepStart, hStep, errmaxSqFallThru, progressedLane && !goodStep, yFinal, hFinal, errmax_sqFinal);
#else
  // Try 5:   2019.05.09   14:40
  cout << " Store Good Values calls - version 5:   2019.05.09  14:40" << std::endl;
  StoreGoodValues(yStepEnd,   hStep,         errmaxSqFallThru, progressedLane &&  goodStep, yFinal, hFinal, errmax_sqFinal);  
  StoreGoodValues(yStepStart, hStepLastGood, errmaxSqLastGood, progressedLane && !goodStep, yFinal, hFinal, errmax_sqFinal);  
#endif
  
  if( 1 ) { 
     cout << " Check LastGood values vs. hstep-last/FallThru" << std::endl;
     ReportRowOfDoubles("hStep =", hStep );
     ReportRowOfDoubles("hStepLastGood = ", hStepLastGood );
     ReportRowOfDoubles("errmaxSqFallThru = ", errmaxSqFallThru );                      
     ReportRowOfDoubles("errmaxSqLastGood = ", errmaxSqLastGood );
  }
  
#ifdef DRIVER_DIAGNOSTICS  
  if (!vecCore::MaskEmpty(UnderflowOccurred)) {
    // Real_v xnewB       = x + hnew;
    int numUnderTot= countMaskTrue<Real_v>( UnderflowOccurred );
    std::cerr << "WARNING > Check after iteration loop: found " << numUnderTot
              << "underflow lanes." << std::endl;
    // ReportConditionLanes(UnderflowOccurred, x, xnewB, hStep, htry);
  }
#endif

#ifdef DRIVER_PRINT_PROGRESS  
  if (partDebug)
    cout << "RollingIntDrv: KeepStepping - Loop done at iter = " << iter
         << " with hdid = " << hdid
         << "  last h = " << hStep
         << "  last good h = " << hStepLastGood
         << " from htry= " << htry
         << std::endl;
#endif
  
  hStep = hFinal;
  // errmax_sq = errmax_sqFinal;

  // Some lanes could have failed - others succeded.  General method needed:
/********
#ifndef CHECK_STRETCH_FACTOR
  // hnext = ComputeNewStepLengthWithinLimits2( errmax_sqFinal, hStep );  
  hnext = ComputeNewStepLengthWithinLimits3( errmax_sqFinal, hStep, isNeeded, errStretchFactor );
#else
  Real_v errStretch = ComputeNewStepLengthWithinLimits2( errmax_sqFinal, Real_v(1.0) );
  hnext = errStretch * hStep;  
#endif
******/
  Real_v  errStretch(1.0);
  Bool_v  isNeeded = ! integrationDone;
  hnext = ComputeNewStepLengthWithinLimits3( errmax_sqFinal, hStep, isNeeded, errStretch );
  
#ifdef CHECK_STRETCH_FACTOR
  // ------------------------------------------
  // The old way (improved) - to cross check
  constexpr const double tinyErr2 = 1e-100;
  Real_v emax2pos          = Max(errmax_sqFinal, Real_v(tinyErr2));
  Real_v powerToUse        = vecCore::Blend( errmax_sqFinal < 1.0, Real_v(GetPowerGrow()), Real_v(GetPowerShrink()) );
  Real_v errStretchOpen    = fSafetyFactor * Exp((0.5 * powerToUse) * Log(emax2pos)); // Was: Log(errmax_sqFinal) )
  // Real_v errStretchOpen     = fSafetyFactor * Pow(emax2pos,(0.5 * GetPowerGrow()));
  // ReportRowOfDoubles( "errStretchOpen (raw)", errStretchOpen);
  const Real_v errStretchMaxed  = Min( errStretchOpen, Real_v(fMaxSteppingIncrease) );
  // ReportRowOfDoubles( "old: errStretch (after min) ", errStretchMaxed);
  // ReportRowOfDoubles( "old: errStretch (after min/ceiling) ", errStretchMaxed);  

  // vecCore::MaskedAssign(errStretchOld, belowMin, Real_v(fMaxSteppingIncrease));
  const Real_v errStretchOld  =   Max( errStretchMaxed, Real_v(fMaxSteppingDecrease) );
  // ReportRowOfDoubles( "old: errStretch (after max) ",   errStretchOld);
  // ReportRowOfDoubles( "old: errStretch (after max/floor) ",   errStretchOld);  

  // Fix the size for steps with zero or very tiny error !!
  // constexpr const double minErr2 = 1e-100;  
  // Bool_v belowMin = errmax_sqFinal <= minErr2;
  // vecCore::MaskedAssign(errStretchOld, belowMin, Real_v(fMaxSteppingIncrease));

  // const double minErr2 = fShrinkThresh * fShrinkThresh ; // Was: 1e-100;
  // vecCore::MaskedAssign(errStretchOld, belowMin, Real_v(fMaxSteppingIncrease));
  
  // ReportRowOfDoubles( "old: errStr (change4tiny) ", errStretchOld);    
  // ReportRowOfDoubles( "old: errStretch", errStretchOld);
  // ------------------- End of Old way -----

  if (!vecCore::MaskEmpty( vecCore::math::Abs(errStretch - errStretchOld) > 0.5e-9 * (errStretch+errStretchOld)))
  {
    cout << "===============---- Error report  START  ----========================================" << std::endl;
    cout << "ERROR> Lanes found with differences in calculated value of 'errStretch'"
         << "       ( factor for stretching step size for 'good' next step." << endl;
    ReportRowOfDoubles("new-old: errStretch", errStretch - errStretchOld);
    ReportRowOfDoubles("old: errStretch", errStretchOld);
    ReportRowOfDoubles("new: errStretch", errStretch);
    cout << "===============---- Error report  END    ----========================================" << std::endl;    
  }
#endif

  // ReportRowOfDoubles( "OGS: h-final", hFinal);

#ifdef DRIVER_PRINT_PROGRESS  
  bool KPSreport = true;
  if ( 1 /*partDebug*/ && KPSreport)
  {
    ReportRowOfDoubles("x: Updated =", x );
// #ifndef CHECK_STRETCH_FACTOR
//     Real_v errStretch = ComputeNewStepLengthWithinLimits2( errmax_sqFinal, Real_v(1.0) );
// #endif
    ReportRowOfDoubles("KPS:  errmax2  ", errmax_sqFinal);
    ReportRowOfSquareRoots("KPS:  errmax   ", errmax_sqFinal);
    ReportRowOfDoubles("KPS:   x (new) ", x);    
    ReportRowOfDoubles("KPS:  h-did    ", hdid);
    // ReportRowOfDoubles("KPS:   x (new) ", x);
    ReportRowOfDoubles("KPS:  h-next   ", hnext);
    ReportRowOfDoubles("KPS: facStretch", errStretch);

    Real_v xFinish = xStart + hdid;
    ReportRowOfDoubles("KPS: x0+hdid ", xFinish);
    ReportRowOfDoubles("KPS: xEnd (arg)", xEnd);
    // ReportRowOfDoubles( "KPS: hFinal", hFinal);
  }
  if (partDebug) {
    cout << " hdid= " << hdid << " and hnext= " << hnext << std::endl
         << " end of  KeepStepping method .......................... " << std::endl;
  }
#endif
  // for(int k=0;k<Nvar ;k++) { y[k] = yFinal[k]; }
  
  return;
} // end of KeepStepping()

// ---------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
inline
   void RollingIntegrationDriver<T_Stepper, Nvar>::
      StoreAllValues( const Real_v                    yWork[],
                      const Real_v                  & hValue,
                      const Real_v                  & errMaxSqValue,
                      Real_v                          yFinal[],
                      Real_v                        & hFinalStore,
                      Real_v                        & errMaxSqStore
         ) const
{
  // std::cout << "StoreGoodValues: Unconditional assignment to output - all together." << std::endl;
   
  for (unsigned int j = 0; j < Nvar; ++j)
     yFinal[j]         = yWork[j];
   
  hFinalStore   = hValue;
  errMaxSqStore = errMaxSqValue;   
}

// ---------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
inline
   void RollingIntegrationDriver<T_Stepper, Nvar>::
      StoreGoodValues( const Real_v                    yWork[],
                       const Real_v                  & hValue,
                       const Real_v                  & errMaxSqValue,
                       const vecCore::Mask_v<Real_v> & storeFlag,
                       Real_v                          yFinal[],
                       Real_v                        & hFinalStore,
                       Real_v                        & errMaxSqStore
         ) const
    // yWork,  hValue,      epsSqVal   represent the  input variables
    // yFinal, hFinalStore, epsSqStore represent the output variables
{
  using std::cout;
  using std::endl;
  using vecCore::MaskedAssign;
  using vecCore::Get;
  using vecCore::Set;

  if (vecCore::MaskFull(storeFlag)) {
     StoreAllValues( yWork, hValue, errMaxSqValue, yFinal, hFinalStore, errMaxSqStore );
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
      // Set( finishedLane, i, true);  //  finishedLane[i] = true;

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
    const FieldTrack yInput[],
    const double     hstep[],
    const double     charge[],
    int              nTracks,
    int       indexArr[], // [VecSize]
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

/******
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
*******/

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
bool RollingIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::InsertNewTrack(
    const FieldTrack yInput[],
    const double hstep[],
    const double charge[],
    const int slot,
    int &trackNextInput,
    bool       succeeded[],      // Output 
    Real_v     y[], // [Nvar]
    Real_v   & hStepLane,
    Real_v   & chargeLane,
    Real_v   & startCurveLength,
    int        indexArr[], // [VecSize]
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
    if (hStepNext > 0) {
      double yScalar[Nvar /* fNoVars*/ ];
      badStepSize[trackNextInput]= false;
      
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
      badStepSize[trackNextInput]= true;
      numBadSize++;
      
      if (hStepNext == 0) {
        std::cerr << "Proposed step is zero; hstep = " << hStepNext << " !" << std::endl;
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
void RollingIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::StoreOutput(const Real_v yEnd[],
                                                                        const Real_v x,
                                                                        int          currIndex,
                                                                        FieldTrack   yOutput[],
                                                                        int          indOut,
                                                                        const double hstep[],
                                                                        bool succeeded[],
                                                                        int nTracks) const
// Called whenever a lane is finished, to store output of integration
// Input arguments
//    - currIndex is the index of finished lane in Vc vector
//    - yEnd[] vector of values of integrands
//    - x      vector of final length (independent parameter value)
//    - hstep  don't store if h<=0
// Stores values into: 
//    - yOutput[indOut]  - the end position and momentum.  ( Note: yOutput is ordinary array. )
//    - yOutput[indOut].fCurveLenght - the final curve length,
//    - succeeded[nTracks]  success flag in array 
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

#ifdef CHECK_ONE_LANE
  const int trackToPrint = IntegrationDriverConstants::GetInstance()->GetTrackToCheck();              
  if( indOut == trackToPrint )
     laneToCheck = -1;
#endif
  
} // End of StoreOutput function

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
void RollingIntegrationDriver<T_Stepper, Nvar>::AccurateAdvance(const FieldTrack yInput[],
                                                                const double     hstep[],
                                                                const double     charge[],
                                                                // double epsilon, // Can be scalar or varying
                                                                FieldTrack yOutput[],
                                                                bool       stillOK[],
                                                                int        nTracks) const
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

  const std::string methodName = "RollingIntgrDriver::AccurateAdvance";
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
  Real_v x, hnext, hdid=0.0, chargeLane, xStartLane, hStepLane;
  Real_v y[Nvar /*ncompSVEC*/ ];     // Working array 1
  Real_v yNext[Nvar /*ncompSVEC*/ ], // Working array 2
      dydx[Nvar /*ncompSVEC*/ ];

  // Real_v ySubStepStart[Nvar /*ncompSVEC*/ ];
  bool *badStepSize = new bool[nTracks];      // Step size is illegal, ie. h <= 0

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
    indexArr[i]       = -1;

  int idNext = InitializeLanes(yInput, hstep, charge, nTracks,                             // Input
                               indexArr, y, hStepLane, chargeLane, xStartLane, numFilled,  // Main output
                               badStepSize, numBadSize);                                   // 'Bad' output
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
    for (unsigned int i = 0; i < VecSize; ++i) // VecSize; ++i)
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
  // ... ELSE just compute this, and update later.
  // Real_v momStart = GetMomentumMag(y);  // True starting momentum
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

  while (
      (   !vecCore::MaskFull(isDoneLane)
          && !vecCore::MaskEmpty((nstp <= fMaxNoSteps)
                                 && (x < x2) && (!isLastStepLane)))
      ||
          idNext < nTracks
     ) {
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
    fpStepper ->RightHandSideInl( y, chargeLane, dydx ); // TODO: change to inline
    //---------****************-----------------------


    ReportRowOfDoubles("hTry(before)", hTry);
    bool lastBatch = ( idNext == nTracks);

    std::cout << " Just before call to KeepStepping:  lastBatch = " << lastBatch << " idNext = " << idNext << " nTracks= "  << nTracks << std::endl;
    
    KeepStepping<Real_v>( y,
                          dydx,        // Returns last value !
                          chargeLane,
                          lastBatch, 
                          x,
                          x2, /* new argument : xEnd */ 
                          hTry,
                          // epsilon,
                          yNext,
                          hdid,
                          // hLast,
                          hnext );  // hStepLane, hDone) ;
    //**********---------------------------------------------------------

#ifdef DRIVER_PRINT_PROGRESS    
    if (partDebug) {
      cout << "##-------------------------------------------------------------------------------" << endl;       
      cout << "### Accurate Advance ---- After return from KeepStepping" << endl;
      ReportRowOfDoubles("charge", chargeLane);
      // ReportRowOfDoubles("h-ask", h);
      ReportRowOfDoubles("hTry (after)", hTry);
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
      // momStart = GetMomentumMag(yStart); // True initial momentum magnitude      
      Real_v momStepStart = GetMomentumMag(y);  //  'Step' start - easier to have
      ReportRowsOfPositionsMomenta("yNext", yNext, Nvar, momStepStart);
      // ReportManyRowsOfDoubles("yNext", yNext, Nvar);
      // Real_v momTrueStart = GetMomentumMag(yStart); // True initial momentum magnitude
      // ReportRowsOfPositionsMomenta("yNext", yNext, Nvar, momTryeStart);  // Check vs true initial mag.
      // cout << "##-------------------------------------------------------------------------------" << endl;
    }
#endif

    // lastStepOK = (hdid == h);
    fNoTotalSteps++;

#ifdef DRIVER_PRINT_PROGRESS      
    bool reportMove = true;
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

    const double epsilonRelMax= FlexIntegrationDriver::GetMaxRelativeEpsilon();
    
    // Note: xStartLane must be positive. ( Ok - expect starting value = 0 )
    Real_v stepThreshold           = vecCore::math::Min(epsilonRelMax * hStepLane, fSmallestFraction * (xStartLane + hdid));
    Bool_v avoidNumerousSmallSteps = hTry < stepThreshold;

    // If it is always true for h<=0 --> lastStep is true, hence the lane will be sent to StoreOutput.

    isLastStepLane = avoidNumerousSmallSteps || isLastStepLane;
    // 'Accumulate' ie use 'OR' - in case lane is already finished or empty

    // x += hdid;  It is already updated - do not add the step again!!

    Real_v xRemains = x2 - x; // (hStepLane - x) + x1;
    // For rest, check the proposed next stepsize

#ifdef DRIVER_PRINT_PROGRESS
    if (partDebug) {
      cout << " hRequest= " << hStepLane << endl;
      // cout << " x-Start = " << xStartLane << endl;
      cout << " hdid    = " << hdid << endl;
      cout << " x (Now) = " << x << endl;
      cout << " x2 -x   = " << xRemains << endl;
    }
#endif

    hnext = Max(hnext, Real_v(GetMinimumStep()));
    // Note: This has potential for instability i.e. if MinStep is 'too long'

    // Ensure that the next step does not overshoot
    hnext = Min(xRemains, hnext);

#ifdef DRIVER_PRINT_PROGRESS
    if ( 1 /*partDebug*/ ) {
       cout << "AccurateAdvance: hnext = " << hnext << " to replace hTry = " << hTry << endl;
       ReportRowOfDoubles("SID/AccAdv: hNext=", hnext);
    }
#endif


    hTry = hnext;
    
    // When stepsize overshoots, decrease it!
    // Must cope with difficult rounding-error issues if hstep << x2

    isLastStepLane = (hTry == 0.0) || isLastStepLane;

    if (partDebug) std::cout << " lastStep : " << isLastStepLane << std::endl;

    nstp += 1; // nstp++;

    succeededLane = (xRemains <= 0.0); // (x>=x2); // If it was a "forced" last step ?

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
      for (unsigned int i = 0; i < VecSize; ++i) {
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
                                         xStartLane, indexArr, nTracks,
                                         badStepSize,  numBadSize );
            
            Set(   isDoneLane, i, !filled);  // isDoneLane[i]   = !filled;
            Set( finishedLane, i, !filled);  // finishedLane[i] = !filled;
            Set( renewedLanes, i,  filled);

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
      } // for ( uint i = 0; i < VecSize; ++i)
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
        MaskedAssign(nstp, renewedLanes, Index_v(0)); // ==> Requires compatible Integer type ...
        // MaskedAssign(  nstp,  renewedLanes, Real_v(0.0) );
        MaskedAssign( x1,   renewedLanes, xStartLane);
        MaskedAssign( x,    renewedLanes, xStartLane);
        MaskedAssign( hTry, renewedLanes, hStepLane);
        // MaskedAssign(  x,  renewedLanes, x1         );       // Done at top of loop - for now
        MaskedAssign( hdid, renewedLanes, Real_v(0.0)); // Maybe not really needed ... but good!

        x2 = x1 + hStepLane; // It does not hurt to do it again for some lanes ...

        // Copy the remaining values ...
        // for (unsigned int i = 0; i < VecSize; ++i) {
        //    y[i] = vecCore::Blend( renewedLanes, yInput[i], yNext[i] );
        //    vecCore::MaskedAssign( y, !renewedLanes, yNext[i] );
        // }

        // vecCore::MaskedAssign( isLastStepLane, renewedLanes, Bool_v(false) );
        isLastStepLane = isLastStepLane && !renewedLanes;

#ifdef DRIVER_PRINT_PROGRESS        
        if (partDebug) {
          cout << " --Inserting New Track - After 'masked' reset of remaining state: " << endl;           
          cout << " *** Vectors changed together after 'loop':" << endl;
          ReportRowOfDoubles("x1", x1);
          ReportRowOfDoubles("x", x);
          ReportRowOfDoubles("hTry", hTry);
          ReportRowOfDoubles("hdid", hdid);
          ReportRowOfDoubles("numStep", nstp);
          // ReportRowOfDoubles( "hDone",   hDone);
          cout << " *** Existing values - recall" << endl;
          ReportRowOfBools<Real_v>("renewedLanes", renewedLanes);
          ReportRowOfDoubles("xStart", xStartLane);          
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
        ReportRowOfDoubles("h(next)", hnext);
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
          for (int i = 0; i < VecSize; ++i)
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
          fpScalarDriver->AccurateAdvance(y_input, hstep[ indexArr[indLastLane] ] - hDone[indLastLane], // epsilonRelMax,
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

  delete[] badStepSize;  // For now - move it to thread member for performance
  
  return;
} // end of AccurateAdvance (flexible / vector )  .....................

//----------------------------------------------------------------------


template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
bool RollingIntegrationDriver<T_Stepper, Nvar>::CheckOutput(Real_v Output[], int lane, int initialSlot,
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
bool RollingIntegrationDriver<T_Stepper, Nvar>::TestInitializeLanes() // int numTracks)
{
  using std::cout;
  using std::cerr;
  using std::endl;
  using vecCore::Get;
  constexpr unsigned int VecSize = vecCore::VectorSize<Real_v>();
  
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

  const int minSize = std::min((int)VecSize, (int)numTracks);
  if (nFilled1 != minSize) { // Can cope with 2, 4 or 8
    cerr << testName << testId << " ERROR> Found less than all lanes filled: Number = " << nFilled1
         << " != expected = " << VecSize << endl;
  }
  assert(nFilled1 == VecSize);

  if (numNext1 != VecSize) { // Can cope with 2, 4 or 8
    cerr << testName << testId << " ERROR> Found next track number = " << numNext1
         << " != expected = " << VecSize << endl;
  }
  assert(numNext1 == VecSize);

  if (numBadSize != 0) {
    cerr << testName << testId << " ERROR> Found non-zero bad size lanes." << endl;
  }
  assert(numBadSize == 0);

  // Calling Checking method
  cout << "      calling Checking method" << endl;
  std::string varNameOutput("yWorkLanes");
  for (unsigned int iSlot = 0; iSlot < VecSize; iSlot++) {
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
  int expectNext2       = predictedLast2[std::min((int)VecSize, 4)];
  if (nextLocation2 != expectNext2) {
    cerr << testName << testId << " ERROR> Found less than the expect id of next location: Number = " << nFilled2
         << "  expected = " << expectNext2 << endl;
  }
  int expectFill2 = std::min((int)VecSize, 3);
  if (nFilled2 != expectFill2) {
    cerr << testName << testId << " ERROR> Found less than the expect number of lanes filled: Number = " << nFilled2
         << "  expected = " << expectFill2 << endl;
  }
  // assert( nFilled2 == std::min( VecSize, 3 ) );
  // int identityNums2[9]= { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
  int predictedNumBad2[9] = {0, 2, 3, 5, 5, 5, 5, 5, 5};
  int expectedNumBad2 =
      predictedNumBad2[VecSize]; // ==? VecSize - nFilled2;
  if (numBadSize != expectedNumBad2)
    cerr << testName << testId << " ERROR> Found " << numBadSize << " bad step-size"
         << " versus " << expectedNumBad2 << " expected for VecSize = " << VecSize << endl;

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
  if (VecSize > 2) {
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
                                                                // double epsilon, // Can be scalar or varying
                                                               FieldTrack yOutput[], int nTracks, bool stillOK[]) const
{
  AccurateAdvance<geant::Double_v>(yInput, hstep, charge,
                                   // epsilon, // Can be scalar or varying
                                   yOutput, stillOK, nTracks);
}

#ifdef EXTEND_SINGLE
template <class T_Stepper, unsigned int Nvar>
void RollingIntegrationDriver<T_Stepper, Nvar>::AccurateAdvance(const FieldTrack & yInput,
                                                                const double       hstep,
                                                                const double       charge, // double epsilon,
                                                                FieldTrack       & yOutput,
                                                                bool               succeeded) const
{
  AccurateAdvance<double>(&yInput, &hstep, &charge,
                          // epsilon, // Can be scalar or varying
                          &yOutput, &succeeded, 1);
}
#endif

#endif /* RollingIntegrationDriver_Def */
