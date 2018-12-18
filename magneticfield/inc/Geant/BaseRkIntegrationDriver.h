//
// Base for Drivers that use explicit Runge-Kutta methods
//
// Class description:
//
// Provides common methods needed for RK drivers
//
// History:
//
// Created: J. Apostolakis,     Nov 2018
//
// Adapted from Simple Integration Driver class
//
// - Contributors: J.Apostolakis                    2018
// ---------------------------------------------------------------

#ifndef BaseRkIntegrationDriver_Def
#define BaseRkIntegrationDriver_Def

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
#include "Geant/FormattedReporter.h"

#include "Geant/VectorTypes.h" //  Defines geant::Double_v
#include "Geant/math_wrappers.h"

// --------------------------------------------------------------
template <class T_Stepper, unsigned int Nvar>
   class BaseRkIntegrationDriver : public FlexIntegrationDriver // , public vecgeom::AlignedBase
{
public:
  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

  BaseRkIntegrationDriver(double     hminimum, // same
                          T_Stepper *pStepper,
                          int        numberOfComponents = 6,
                          int        statsVerbosity = 1);

  virtual ~BaseRkIntegrationDriver();

  // Setting parameters ( few now )
  void SetMaxNoSteps(int val) { fMaxNoSteps = val; }
  unsigned long IncrementStepperCalls() const { return ++fStepperCalls; } // mutable ..
  unsigned long GetNumberOfStepperCalls()    { return fStepperCalls; }

  template <typename Real_v>
  Real_v ComputeNewStepSize( Real_v errMaxNorm,   // max error  (normalised)
                             Real_v hStepCurrent);  // current step size

  template <typename Real_v>  
  Real_v ComputeNewStepSize_WithinLimits(
    Real_v errMaxNorm,   // max error  (normalised)
    Real_v hStepCurrent); // current step size

  /***
  inline const T_Stepper *GetStepper() const { return fpStepper; }
  inline       T_Stepper *GetStepper()       { return fpStepper; }
   ***/
  
  double GetPowerShrink()  const { return fPowerShrink; }
  double GetPowerGrow()    const { return fPowerGrow; }
  int    GetVerboseLevel() const { return fVerboseLevel; }
  int    GetErrcon()       const { return fErrcon; }  
  int    GetStatisticsVerboseLevel() const { return fStatisticsVerboseLevel; }
  double GetMinimumStep()  const { return fMinimumStep; }
  unsigned int GetStepperOrder() const { return fStepperOrder; }
  // Accessors.

  void SetVerboseLevel( int val ) { fVerboseLevel= val; }
  void SetMinimumStep(double dval) { fMinimumStep = dval; }
  // Modifiers
  
  // Compute dependent parameters
  inline void ComputeAndSetErrcon();

  void ReportInvalidStepInputs(double hStepArr[], int nTracks);
  void CreateInput(FieldTrack yInput[], int nTracks);

  // int  BookStatistics();
  // int  ReportStatistics() const;
  
  // Checks
  void CheckParameters(); // Sanity check of
  template <typename Real_v>  
  bool CheckOutput(Real_v Output[], int lane, int initialSlot, std::string testName, std::string varName);
      
protected:  //  Copy constructor and no assignment operator.

  BaseRkIntegrationDriver(const BaseRkIntegrationDriver &);
  // Copy constructor used to create Clone method

private:
  BaseRkIntegrationDriver &operator=(const BaseRkIntegrationDriver &) = delete;

protected:
  // ---------------------------------------------------------------
  // DEPENDENT Objects
  // T_Stepper *fpStepper;
  
  // ---------------------------------------------------------------
  //  INVARIANTS

  double fMinimumStep; // same
  // Minimum Step allowed in a Step (in absolute units)
  double fSmallestFraction = 1.0e-7; // Expected value: larger than 1e-12 to 5e-15;
  // Smallest fraction of (existing) curve length - in relative units
  //  below this fraction the current step will be the last

  // const int  fNoIntegrationVariables;  // Number of Variables in integration
  const int fMinNoVars; // Minimum number for TemplateFieldTrack<Real_v>
  const int fNoVars;    // Full number of variable

  unsigned long fMaxNoSteps;
  static constexpr int fMaxStepBase = 250;

  // static constexpr double fSafetyFactor= 0.9; // -> Fails to compile on clang 9.1 2017.12.05
  const double fSafetyFactor = 0.9; //     OK ...
  const unsigned int fStepperOrder = 0; //  the Integrator order of the (RK or other) stepper
  const double fPowerShrink;        //  exponent for shrinking
  const double fPowerGrow;          //  exponent for growth
  /*const*/ double fErrcon;
  // Parameters used to grow and shrink trial stepsize.

  // double fSurfaceTolerance = 1.e-6;

  //  Stepsize can increase by no more than 5.0
  //           and decrease by no more than x10. = 0.1
  static constexpr double fMaxSteppingIncrease = 5.0;
  static constexpr double fMaxSteppingDecrease = 0.1;
  // Maximum stepsize increase/decrease factors.

  int fStatisticsVerboseLevel         = 0;
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

  mutable unsigned long fStepperCalls = 0UL;
  mutable unsigned long fNoTotalSteps = 0, fNoBadSteps = 0, fNoSmallSteps = 0, fNoInitialSmallSteps = 0;

#ifdef STATISTICS_DEV
  mutable Real_v fDyerr_max(0.0), fDyerr_mx2(0.0);
  mutable Real_v fDyerrPos_smTot(0.0), fDyerrPos_lgTot(0.0), fDyerrVel_lgTot(0.0);
  mutable Real_v fSumH_sm, fSumH_lg;
// Step Statistics
#endif

  int fVerboseLevel; // Verbosity level for printing (debug, ..)
                     // Could be varied during tracking - to help identify issues

}; // End of class definition -- BaseRkIntegrationDriver

/*****
template <class Real_v, class T_Stepper, unsigned int Nvar>
constexpr double BaseRkIntegrationDriver<Real_v, T_Stepper, Nvar>::fMaxSteppingIncrease;

template <class Real_v, class T_Stepper, unsigned int Nvar>
constexpr double BaseRkIntegrationDriver<Real_v, T_Stepper, Nvar>::fMaxSteppingDecrease;
 *****/

// ---------------------------------------------------------
//  Constructor
//
template <class T_Stepper, unsigned int Nvar>
BaseRkIntegrationDriver<T_Stepper, Nvar>::BaseRkIntegrationDriver(double     hminimum,
                                                                  T_Stepper *pStepper,
                                                                  int        numComponents,
                                                                  int        statisticsVerbose)
    :
      fMinimumStep(hminimum),
      // -- fSmallestFraction( 1.0e-12 ),
      // fNoIntegrationVariables(numComponents),  // ==> Nvar
      fMinNoVars(6), fNoVars(Nvar), // Was fNoVars(std::max((int)Nvar, std::max((int)fMinNoVars, (int)numComponents))),
      fStepperOrder( pStepper->GetIntegratorOrder() ),
      fPowerShrink(-1.0 / fStepperOrder),       //  exponent for shrinking
      fPowerGrow(-1.0 / (1.0 + fStepperOrder)), //  exponent for growth
      // - fErrcon(0.0),
      fStatisticsVerboseLevel(statisticsVerbose),
      fVerboseLevel(0)
      // fNoInitialSmallSteps(0),
{
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
  assert(pStepper != nullptr);
  assert(Nvar <= (unsigned int)numComponents); // Ensure that arrays are large enough for Integr.

  // fpStepper = pStepper;

  ComputeAndSetErrcon();
  SetMaxNoSteps( fMaxStepBase / pStepper->GetIntegratorOrder() );

  ComputeAndSetErrcon();

  CheckParameters();

#ifdef GUDEBUG_FIELD
  fVerboseLevel = 2;
#endif

  if (fVerboseLevel) {
    std::cout << "SiD:ctor> Stepper Order= " << pStepper->GetIntegratorOrder() << " > Powers used: "
              << " shrink = " << fPowerShrink << "  grow = " << fPowerGrow << std::endl;
  }
  if ((fVerboseLevel > 0) || (fStatisticsVerboseLevel > 1)) {
    std::cout << "BaseRkIntegrationDriver constructor called.";
       //     << "invE_nS, QuickAdv-2sqrt "
       // << "with Statistics " << fStatsStatus << std::endl;
    // ( fStatsEnabled ? "enabled" : "disabled" )
  }
}

// ---------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
inline void BaseRkIntegrationDriver<T_Stepper, Nvar>
    // void BaseRkIntegrationDriver<Real_v, T_Stepper, Nvar>
    ::ComputeAndSetErrcon()
{
  fErrcon = Math::Pow(fMaxSteppingIncrease / fSafetyFactor, 1.0 / fPowerGrow);
  // return fErrcon;
}

// ---------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
inline void BaseRkIntegrationDriver<T_Stepper, Nvar>::CheckParameters()
{
  constexpr double perMillion = 1.0e-6;
  using std::cerr;
  using std::endl;

  double checkPowerShrink = -1.0 / fStepperOrder;

  double diffShrink = fPowerShrink - checkPowerShrink;
  if (std::fabs(diffShrink) // checkPowerShrink - fPowerShrink)
      >= perMillion * std::fabs(fPowerShrink)) {
    cerr << "BaseRkIntegrationDriver: ERROR in fPowerShrink" << std::endl;
    cerr << "    calculated = " << checkPowerShrink << "    pre-computed = " << fPowerShrink << "  diff= " << diffShrink
         << "  tolerance = " << perMillion * std::fabs(fPowerShrink) << endl;
    cerr << "  Order of integrator = " << fStepperOrder << endl;
    exit(1);
  }
  assert(std::fabs(checkPowerShrink - fPowerShrink) < perMillion * std::fabs(fPowerShrink));

  double checkPowerGrow = -1.0 / (1.0 + fStepperOrder);
  assert(std::fabs(checkPowerGrow - fPowerGrow) < perMillion * std::fabs(fPowerGrow));

  if (std::fabs(checkPowerGrow - fPowerGrow) >= perMillion * std::fabs(fPowerGrow)) {
    std::cerr << "BaseRkIntegrationDriver: ERROR in fPowerGrow" << std::endl;
    exit(1);
  }

  if (fVerboseLevel)
    std::cout << "BaseRkIntegrationDriver::CheckParameters > Powers used: " << std::endl
              << "  shrink = " << fPowerShrink << "  grow = " << fPowerGrow << std::endl;
}

// ---------------------------------------------------------

//  Copy Constructor - used by Clone
template <class T_Stepper, unsigned int Nvar>
BaseRkIntegrationDriver<T_Stepper, Nvar>::BaseRkIntegrationDriver(
    const BaseRkIntegrationDriver</*Real_v,*/ T_Stepper, Nvar> &right)
    : fMinimumStep(right.fMinimumStep), fSmallestFraction(right.fSmallestFraction),
      // fNoIntegrationVariables( right.fNoIntegrationVariables ),
      fMinNoVars(right.fMinNoVars), fNoVars(std::max((int)Nvar, fMinNoVars)), fPowerShrink(right.fPowerShrink),
      fPowerGrow(right.fPowerGrow), fErrcon(right.fErrcon),
      // fSurfaceTolerance( right.fSurfaceTolerance ),
      fStatisticsVerboseLevel(right.fStatisticsVerboseLevel),
      /* fDyerr_max(0.0), fDyerr_mx2(0.0),
         fDyerrPos_smTot(0.0), fDyerrPos_lgTot(0.0), fDyerrVel_lgTot(0.0),
         fSumH_sm(0.0), fSumH_lg(0.0), */
      fVerboseLevel(right.fVerboseLevel)
{
  // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
  // is required. For proper time of flight and spin,  fMinNoVars must be 12
  // const T_Stepper *protStepper = right.GetStepper();
  // fpStepper                    = protStepper->Clone();

  ComputeAndSetErrcon();
  fMaxNoSteps = fMaxStepBase / fStepperOrder; // fpStepper->GetIntegratorOrder();

  if ((fVerboseLevel > 0) || (fStatisticsVerboseLevel > 1)) {
    std::cout << "BaseRkIntegrationDriver copy constructor called."
           // << "with Statistics flag = " << fStatsStatus
              << std::endl;
  }
}

// ---------------------------------------------------------

//  Destructor
template <class T_Stepper, unsigned int Nvar>
BaseRkIntegrationDriver<T_Stepper, Nvar>::~BaseRkIntegrationDriver()
{
  if (fStatisticsVerboseLevel > 1) {
    // PrintStatisticsReport();
  }

  // delete fpScalarDriver;
  // delete fpScalarStepper;
  // delete fpStepper;
}

// --------------------------------------------------------------------------

//  This method computes new step sizes - but does not limit changes to
//  within  certain factors
//
template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
Real_v BaseRkIntegrationDriver</*Real_v,*/ T_Stepper, Nvar>::
   ComputeNewStepSize(
    Real_v errMaxNorm,   // max error  (normalised)
    Real_v hStepCurrent) // current step size
{
  using Bool_v = vecCore::Mask_v<Real_v>;

  Bool_v goodStep = (errMaxNorm <= 1.0);
  Real_v powerUse = vecCore::Blend(goodStep, fPowerShrink, fPowerGrow);
  Real_v stretch  = fSafetyFactor * vecgeom::Pow(errMaxNorm, powerUse);
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
Real_v BaseRkIntegrationDriver<T_Stepper, Nvar>::
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

  /*
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
    }*/

  return hNew;
}

// ---------------------------------------------------------------------------
template <class T_Stepper, unsigned int Nvar>
void BaseRkIntegrationDriver<T_Stepper, Nvar>::ReportInvalidStepInputs(double hStepArr[], int nTracks)
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

// ####################  Testing method(s) ####################

#include <cassert>

template <class T_Stepper, unsigned int Nvar>
void BaseRkIntegrationDriver<T_Stepper, Nvar>::CreateInput(FieldTrack yInput[], int nTracks)
{
  constexpr int Ncomp = FieldTrack::NumCompFT;
  double PositionMomArr[Ncomp];

  assert(Nvar <= Ncomp); // Ensure that arrays are large enough for Integr.

  double length;
  for (int itr = 0; itr < nTracks; itr++) {
    for (unsigned int j = 0; j < Nvar; j++) {
      PositionMomArr[j] = 10. * itr + j;
    }
    length = 0.1 * itr;
    // yInput[j] = FieldTrack( PositionMomArr, length );
    yInput[itr].LoadFromArray(PositionMomArr, Nvar);
    yInput[itr].SetCurveLength(length);
  }
}

// ------------------------------------------------------------

template <class T_Stepper, unsigned int Nvar>
template <class Real_v>
bool BaseRkIntegrationDriver<T_Stepper, Nvar>::
   CheckOutput(Real_v Output[], int lane, int initialSlot, std::string testName, std::string varName)
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

// template <class T_Stepper, unsigned int Nvar>
// int BaseRkIntegrationDriver<T_Stepper, Nvar>::
//    ReportStatistics() const
// {
// }

// ------------------------------------------------------------
// #include <TH1.h>

// template <class T_Stepper, unsigned int Nvar>
// int BaseRkIntegrationDriver<T_Stepper, Nvar>::BookStatistics() 
// {
// }

#endif /* BaseRkIntegrationDriver_Def */
