//
//
// Derived from G4MagInt_Drv class of Geant4 (in G4MagIntegratorDriver.hh/cc)
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and insures that 
// the error is within acceptable bounds.
//
// History:
// - Created. J.Apostolakis.
// --------------------------------------------------------------------

#ifndef ScalarIntegrationDriver_Def
#define ScalarIntegrationDriver_Def

#include "ScalarFieldTrack.h"

// class GUVIntegrationStepper;
#include "GUVIntegrationStepper.h"

class ScalarIntegrationDriver
{
   public:  // with description
     ScalarIntegrationDriver( double                 hminimum, //same
                          GUVIntegrationStepper *pStepper,
                          int                    numberOfComponents=6,
                          int                    statisticsVerbosity=1);
     ScalarIntegrationDriver( const ScalarIntegrationDriver& );
       // Copy constructor used to create Clone method
     ~ScalarIntegrationDriver();

     // Core methods
     bool  AccurateAdvance( const ScalarFieldTrack& y_current,
                                        double  hstep,
                                        double  eps, //same             // Requested y_err/hstep
                                  ScalarFieldTrack& yOutput,                            
                                        double  hinitial=0.0);  // take it out 
       // Above drivers for integrator (Runge-Kutta) with stepsize control. 
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps. 
       // On output ystart is replaced by value at end of interval. 
       // The concept is similar to the odeint routine from NRC p.721-722.

     bool  QuickAdvance(      ScalarFieldTrack& y_posvel,        // INOUT
                          const double      dydx[],  
                                double      hstep,           // IN
#ifdef USE_DCHORD
                                double&     dchord_step,
#endif
                                double&     dyerr_pos_sq,
                                double&     dyerr_mom_rel_sq ) ;
       // New QuickAdvance that also just tries one Step
       //    (so also does not ensure accuracy)
       //    but does return the errors in  position and
       //        momentum (normalised: Delta_Integration(p^2)/(p^2) )

     ScalarIntegrationDriver* Clone() const;
       // Create an independent copy of the current object -- including independent 'owned' objects
       // 
       // Question:  If the current object and all sub-objects are const, can it return 'this' ?
     
     GUVEquationOfMotion* GetEquationOfMotion() { return fpStepper->GetEquationOfMotion(); }
     const GUVEquationOfMotion* GetEquationOfMotion() const { return fpStepper->GetEquationOfMotion(); }

     // Auxiliary methods
     inline double GetHmin()        const { return fMinimumStep; } 
     inline double GetSafety()      const { return fSafetyFactor; }
     inline double GetPowerShrink() const { return fPowerShrink; }
     inline double GetPowerGrow()   const { return fPowerGrow; } 
     inline double GetErrcon()      const { return fErrcon; }
     
     inline void   GetDerivatives( const ScalarFieldTrack &y_curr,     // const, INput
                                        double    charge,
                                        double    dydx[]   );  //       OUTput
        // Accessors.

     inline void RenewStepperAndAdjust(GUVIntegrationStepper *Stepper);
        // Sets a new stepper 'Stepper' for this driver. Then it calls
        // ReSetParameters to reset its parameters accordingly.

     inline void ReSetParameters(double new_safety= 0.9 );
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

     inline const GUVIntegrationStepper* GetStepper() const;
     inline GUVIntegrationStepper* GetStepper();

     void  OneGoodStep(       double  ystart[], // Like old RKF45step()
                              double  charge,
                        const double  dydx[],
                              double& x,
                              double htry,
                              double  eps,      //  memb variables ?
                              double& hdid,
                              double& hnext ) ;
        // This takes one Step that is as large as possible while 
        // satisfying the accuracy criterion of:
        // yerr < eps * |y_end-y_start|

     double ComputeNewStepSize( double  errMaxNorm,    // normalised error
                                  double  hstepCurrent); // current step size
        // Taking the last step's normalised error, calculate
        // a step size for the next step.
        // Do not limit the next step's size within a factor of the
        // current one.

     double ComputeNewStepSize_WithinLimits(
                          double  errMaxNorm,    // normalised error
                          double  hstepCurrent); // current step size
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
     void WarnSmallStepSize( double hnext, double hstep, 
                             double h,     double xDone,
                             int noSteps);
     void WarnTooManySteps( double x1start, double x2end, double xCurrent);
     void WarnEndPointTooFar (double  endPointDist, 
                              double  hStepSize , 
                              double  epsilonRelative,
                              int     debugFlag);
        //  Issue warnings for undesirable situations

     void PrintStatus(  const double*      StartArr,
                              double       xstart,
                        const double*      CurrentArr, 
                              double       xcurrent, 
                              double       requestStep, 
                              int          subStepNo );
     void PrintStatus(  const ScalarFieldTrack&  StartFT,
                        const ScalarFieldTrack&  CurrentFT, 
                              double       requestStep, 
                              int          subStepNo );
     void PrintStat_Aux( const ScalarFieldTrack& aScalarFieldTrack,
                               double      requestStep, 
                               double      actualStep,
                               int         subStepNo,
                               double      subStepSize,
                               double      dotVelocities );       
       //  Verbose output for debugging

     void PrintStatisticsReport() ;
       //  Report on the number of steps, maximum errors etc.

#ifdef QUICK_ADV_TWO
     bool QuickAdvance(      double     yarrin[],     // In
                         const double     dydx[],  
                               double     hstep,        
                               double     yarrout[],    // Out
                               double&    dchord_step,  // Out
                               double&    dyerr );      // in length
#endif

   private:

     ScalarIntegrationDriver& operator=(const ScalarIntegrationDriver&);
        // Private copy constructor and assignment operator.

   private:

     // ---------------------------------------------------------------
     // DEPENDENT Objects
     GUVIntegrationStepper *fpStepper;
     
     // ---------------------------------------------------------------
     //  INVARIANTS 

     double  fMinimumStep;
        // Minimum Step allowed in a Step (in absolute units)
     double  fSmallestFraction;      //   Expected range 1e-12 to 5e-15;  
        // Smallest fraction of (existing) curve length - in relative units
        //  below this fraction the current step will be the last 

     const int  fNoIntegrationVariables;  // Number of Variables in integration
     const int  fMinNoVars;               // Minimum number for ScalarFieldTrack
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
     int  fNoTotalSteps, fNoBadSteps, fNoSmallSteps, fNoInitialSmallSteps;
     double fDyerrPosMaxSq, fDyerrDirMaxSq;
     double fDyerrPos_smTot, fDyerrPos_lgTot, fDyerrVel_lgTot;
     double fSumH_sm, fSumH_lg;
        // Step Statistics

     int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
        // Could be varied during tracking - to help identify issues
     int fStepperCalls=0.;     

};

// #include "GUIntegratorDriver.icc"




inline
double ScalarIntegrationDriver::ComputeAndSetErrcon()
{
      fErrcon = std::pow(fMaxSteppingIncrease/fSafetyFactor,1.0/fPowerGrow);
      return fErrcon;
} 

inline
void ScalarIntegrationDriver::ReSetParameters(double new_safety)
{
      fSafetyFactor = new_safety;
      fPowerShrink  = -1.0 / fpStepper->IntegratorOrder();
      fPowerGrow    = -1.0 / (1.0 + fpStepper->IntegratorOrder());
      ComputeAndSetErrcon();
}

inline
void ScalarIntegrationDriver::SetSafety(double val)
{ 
      fSafetyFactor=val;
      ComputeAndSetErrcon();
}

inline
void ScalarIntegrationDriver::SetPowerGrow(double  val)
{ 
      fPowerGrow=val;
      ComputeAndSetErrcon(); 
}

inline
void ScalarIntegrationDriver::SetErrcon(double val)
{ 
      fErrcon=val;
}

inline
void ScalarIntegrationDriver::RenewStepperAndAdjust(GUVIntegrationStepper *pStepper)
{  
      fpStepper = pStepper; 
      ReSetParameters();
}

inline
const GUVIntegrationStepper* ScalarIntegrationDriver::GetStepper() const
{
  return fpStepper;
}

inline
GUVIntegrationStepper* ScalarIntegrationDriver::GetStepper() 
{
  return fpStepper;
}

inline
int ScalarIntegrationDriver::GetMaxNoSteps() const
{
  return fMaxNoSteps;
}

inline
void ScalarIntegrationDriver::SetMaxNoSteps(int val)
{
  fMaxNoSteps= val;
}

inline
void ScalarIntegrationDriver::GetDerivatives(const ScalarFieldTrack &y_curr, // const, INput
                                               double       charge,
                                               double       dydx[])  // OUTput
{ 
  double  tmpValArr[ScalarFieldTrack::ncompSVEC];
  y_curr.DumpToArray( tmpValArr  );
  fpStepper -> RightHandSideVIS( tmpValArr , charge, dydx );
}
#endif /* ScalarIntegrationDriver_Def */
