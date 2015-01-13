//
//
// Derived from GVIntegratorDriver class of Geant4
//
// Class description:
//
// Provides a driver that talks to the Integrator Stepper, and insures that 
// the error is within acceptable bounds.

// History:
// - Created. J.Apostolakis.
// --------------------------------------------------------------------

#ifndef GUIntegrationDriver_Def
#define GUIntegrationDriver_Def

// #include "G4Types.hh"
#include "GUFieldTrack.h"

// class GUVIntegrationStepper;
#include "GUVIntegrationStepper.h"

class GUIntegrationDriver
{
   public:  // with description
     GUIntegrationDriver( double            hminimum, 
                      GUVIntegrationStepper *pItsStepper,
                      int                   numberOfComponents=6,
                      int                   statisticsVerbosity=1);
     ~GUIntegrationDriver();

     // Core methods
     bool  AccurateAdvance( GUFieldTrack& y_current,
                                    double  hstep,
                                    double  eps,            // Requested y_err/hstep
                                    double  hinitial=0.0);  // Suggested 1st interval
       // Above drivers for integrator (Runge-Kutta) with stepsize control. 
       // Integrates ODE starting values y_current
       // from current s (s=s0) to s=s0+h with accuracy eps. 
       // On output ystart is replaced by value at end of interval. 
       // The concept is similar to the odeint routine from NRC p.721-722.

     bool  QuickAdvance(      GUFieldTrack& y_posvel,        // INOUT
                          const double      dydx[],  
                                double      hstep,           // IN
                                double&     dchord_step,
                                double&     dyerr_pos_sq,
                                double&     dyerr_mom_rel_sq ) ;
       // New QuickAdvance that also just tries one Step
       //    (so also does not ensure accuracy)
       //    but does return the errors in  position and
       //        momentum (normalised: Delta_Integration(p^2)/(p^2) )


     // Auxiliary methods
     inline double GetHmin() const;
     inline double Hmin() const;     // Obsolete
     inline double GetSafety() const;
     inline double GetPowerShrink() const;
     inline double GetPowerGrow() const;
     inline double GetErrcon() const;
     inline void   GetDerivatives( const GUFieldTrack &y_curr,     // const, INput
                                       double    dydx[]   );  //       OUTput
        // Accessors.

     inline void RenewStepperAndAdjust(GUVIntegrationStepper *pItsStepper);
        // Sets a new stepper pItsStepper for this driver. Then it calls
        // ReSetParameters to reset its parameters accordingly.

     inline void ReSetParameters(double new_safety= 0.9 );
        //  i) sets the exponents (fPowerGrow & fPowerShrink), 
        //     using the current Stepper's order, 
        // ii) sets the safety
        // ii) calculates "errcon" according to the above values.

     inline void SetSafety(double valS);
     inline void SetPowerShrink(double valPs);
     inline void SetPowerGrow (double valPg);
     inline void SetErrcon(double valEc);
        // When setting safety or fPowerGrow, errcon will be set to a 
        // compatible value.

     inline double ComputeAndSetErrcon();

     inline const GUVIntegrationStepper* GetStepper() const;
     inline GUVIntegrationStepper* GetStepper();

     void  OneGoodStep(       double  ystart[], // Like old RKF45step()
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
     inline void     SetMaxNoSteps( int val); 
        //  Modify and Get the Maximum number of Steps that can be
        //   taken for the integration of a single segment -
        //   (ie a single call to AccurateAdvance).

   public:  // without description

     inline void SetHmin(double newval);
     inline void SetVerboseLevel(int newLevel); 
     inline double GetVerboseLevel() const;

     inline double GetSmallestFraction() const; 
     void     SetSmallestFraction( double val ); 

   protected:  // without description
     void WarnSmallStepSize( double hnext, double hstep, 
                             double h,     double xDone,
                             int noSteps);
     void WarnTooManyStep( double x1start, double x2end, double xCurrent);
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
     void PrintStatus(  const GUFieldTrack&  StartFT,
                        const GUFieldTrack&  CurrentFT, 
                              double       requestStep, 
                              int          subStepNo );
     void PrintStat_Aux( const GUFieldTrack& aGUFieldTrack,
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

     GUIntegrationDriver(const GUIntegrationDriver&);
     GUIntegrationDriver& operator=(const GUIntegrationDriver&);
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
     static const int  fMaxStepBase;  

     double fSafetyFactor;
     double fPowerShrink;   //  exponent for shrinking
     double fPowerGrow;    //  exponent for growth
     double errcon;
        // Parameters used to grow and shrink trial stepsize.

     double fSurfaceTolerance; 

     static const double max_stepping_increase;
     static const double max_stepping_decrease;
        // Maximum stepsize increase/decrease factors.

     int    fStatisticsVerboseLevel;

     // ---------------------------------------------------------------
     // DEPENDENT Objects
     GUVIntegrationStepper *pIntStepper;

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

// #include "GVIntegratorDriver.icc"

inline
double GUIntegrationDriver::GetHmin() const
{
      return fMinimumStep;
} 

inline
double GUIntegrationDriver::Hmin() const
{
      return fMinimumStep;
}

inline
double GUIntegrationDriver::GetSafety() const
{
      return fSafetyFactor;
}

inline
double GUIntegrationDriver::GetPowerShrink() const
{
      return fPowerShrink;
} 

inline
double GUIntegrationDriver::GetPowerGrow() const
{
      return fPowerGrow;
}
 
inline
double GUIntegrationDriver::GetErrcon() const
{
      return errcon;
}

inline
void GUIntegrationDriver::SetHmin(double newval)
{
      fMinimumStep = newval;
} 

inline
double GUIntegrationDriver::ComputeAndSetErrcon()
{
      errcon = std::pow(max_stepping_increase/fSafetyFactor,1.0/fPowerGrow);
      return errcon;
} 

inline
void GUIntegrationDriver::ReSetParameters(double new_safety)
{
      fSafetyFactor = new_safety;
      fPowerShrink  = -1.0 / pIntStepper->IntegratorOrder();
      fPowerGrow    = -1.0 / (1.0 + pIntStepper->IntegratorOrder());
      ComputeAndSetErrcon();
}

inline
void GUIntegrationDriver::SetSafety(double val)
{ 
      fSafetyFactor=val;
      ComputeAndSetErrcon();
}

inline
void GUIntegrationDriver::SetPowerGrow(double  val)
{ 
      fPowerGrow=val;
      ComputeAndSetErrcon(); 
}

inline
void GUIntegrationDriver::SetErrcon(double val)
{ 
      errcon=val;
}

inline
void GUIntegrationDriver::RenewStepperAndAdjust(GUVIntegrationStepper *pItsStepper)
{  
      pIntStepper = pItsStepper; 
      ReSetParameters();
}

inline
const GUVIntegrationStepper* GUIntegrationDriver::GetStepper() const
{
  return pIntStepper;
}

inline
GUVIntegrationStepper* GUIntegrationDriver::GetStepper() 
{
  return pIntStepper;
}

inline
int GUIntegrationDriver::GetMaxNoSteps() const
{
  return fMaxNoSteps;
}

inline
void GUIntegrationDriver::SetMaxNoSteps(int val)
{
  fMaxNoSteps= val;
}

inline
void GUIntegrationDriver::GetDerivatives(const GUFieldTrack &y_curr, // const, INput
                                           double     dydx[])  // OUTput
{ 
  double  tmpValArr[GUFieldTrack::ncompSVEC];
  y_curr.DumpToArray( tmpValArr  );
  pIntStepper -> RightHandSide( tmpValArr , dydx );
}

inline
double GUIntegrationDriver::GetVerboseLevel() const
{
      return fVerboseLevel;
} 

inline 
void GUIntegrationDriver::SetVerboseLevel(int newLevel)
{
      fVerboseLevel= newLevel;
}

inline
double GUIntegrationDriver::GetSmallestFraction() const
{
      return fSmallestFraction; 
} 


#endif /* GUIntegrationDriver_Def */
