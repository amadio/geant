//
// Runge-Kutta Stepper using Cash Karp's RK tableau
//
// Adapted from 'GUTCashKarpRKF45' by Qieshen Xie, GSoC 2014
//         (derived from G4CashKarpRKF45)
//
// First version:  John Apostolakis,  4 Nov 2015
//
#ifndef TCASHKARPRKF45_H
#define TCASHKARPRKF45_H

#include "GULineSection.h"
#include "GUVIntegrationStepper.h"

#define INLINERHS 1

#ifdef INLINERHS
#define REALLY_INLINE   inline __attribute__((always_inline)) 
#else
#define REALLY_INLINE   inline
#endif

template
<class T_Equation, int Nvar>
class GUTCashKarpRKF45 : public GUVIntegrationStepper
{
  public:  
    static constexpr unsigned int sOrderMethod= 4;
    static constexpr unsigned int sNstore = (GUIntegrationNms::NumVarBase > Nvar) ? GUIntegrationNms::NumVarBase : Nvar;
                        // MaxConst( GUIntegrationNms::NumVarBase,  Nvar);
                        // std::max( GUIntegrationNms::NumVarBase,  Nvar);
    // static const IntegratorCorrection = 1./((1<<4)-1); 
    inline double IntegratorCorrection() { return 1./((1<<4)-1); }
   
  public:
    inline 
    GUTCashKarpRKF45(T_Equation *EqRhs,
                   unsigned int noIntegrationVariables=6,
                   unsigned int numStateVariables=0,
                   bool primary=true);

    inline 
    virtual ~GUTCashKarpRKF45();

    REALLY_INLINE
    void RightHandSideInl(double y[], double dydx[]) 
    {fEquation_Rhs->T_Equation::RightHandSide(y, dydx);}

    REALLY_INLINE void StepWithErrorEstimate(const double* __restrict__ yInput,
                                             const double* __restrict__ dydx,
                                             double Step,
                                             double* __restrict__ yOut, // [],
                                             double* __restrict__ yErr);

    double  DistChord()   const;  

    // inline int IntegratorOrder() const { return 4; }

    private:
        GUTCashKarpRKF45(const GUTCashKarpRKF45&) = delete;
        GUTCashKarpRKF45& operator=(const GUTCashKarpRKF45&) = delete;
        //private copy constructor and assignment operator.

    private:
        // State -- intermediate values used during RK step
        // -----        
        double ak2[sNstore];
        double ak3[sNstore];
        double ak4[sNstore];
        double ak5[sNstore];
        double ak6[sNstore];
        double ak7[sNstore];
        double yTemp[sNstore];
        double yIn[sNstore];
        // scratch space

        // State -- values used for subsequent call to DistChord
        // -----
        double  fLastStepLength;
        double* fLastInitialVector;
        double* fLastFinalVector;
        double* fLastDyDx; 
        double* fMidVector;
        double* fMidError;
        // for DistChord calculations

        // Invariants during integration
        // ----------
        GUTCashKarpRKF45* fAuxStepper; 
        T_Equation* fEquation_Rhs;
};


template <class T_Equation, int Nvar>
inline
GUTCashKarpRKF45<T_Equation,Nvar>::
   GUTCashKarpRKF45(T_Equation *EqRhs,
                     unsigned int noIntegrationVariables,
                     unsigned int numStateVariables,
                bool primary)
   : GUVIntegrationStepper(
            dynamic_cast<GUVEquationOfMotion*>(EqRhs),
            sOrderMethod, 
            noIntegrationVariables,   // Redundant - to be DELETED
            ((numStateVariables>0) ? numStateVariables : sNstore) ),
     fLastStepLength(0.), fAuxStepper(0),
     fEquation_Rhs(EqRhs)
{
   const int numberOfVariables = noIntegrationVariables;

   assert( dynamic_cast<GUVEquationOfMotion*>(EqRhs) != 0 );  
   assert( (numStateVariables == 0) || (numStateVariables >= Nvar) );
      
   fLastInitialVector = new double[sNstore] ;
   fLastFinalVector = new double[sNstore] ;
   fLastDyDx = new double[sNstore];
   
   fMidVector = new double[sNstore];
   fMidError =  new double[sNstore];

   std::cerr << " GUTCashKarpRKF45 - constructed class. " << std::endl
             << " Nvar = " << Nvar << " Nstore= " << sNstore 
             << " Primary = " << primary << std::endl;
   if( primary )
   { 
      fAuxStepper = new GUTCashKarpRKF45(EqRhs, numberOfVariables, numStateVariables, !primary);
   }
}

template <class T_Equation, int Nvar>
// inline
REALLY_INLINE
GUTCashKarpRKF45<T_Equation,Nvar>::~GUTCashKarpRKF45()
{
   delete[] fLastInitialVector;
   delete[] fLastFinalVector;
   delete[] fLastDyDx;
   delete[] fMidVector;
   delete[] fMidError;
   
   delete fAuxStepper;
}

template <class T_Equation, int Nvar>
inline void
GUTCashKarpRKF45<T_Equation,Nvar>::
   StepWithErrorEstimate(const double* __restrict__ yInput, // [],    
                         const double* __restrict__ dydx, // [],
                         double Step,
                         double* __restrict__ yOut, // [],
                         double* __restrict__ yErr) // [])
{
    // const int nvar = 6 ;
    // const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
    int i;

    const double  b21 = 0.2 ,
          b31 = 3.0/40.0 , b32 = 9.0/40.0 ,
          b41 = 0.3 , b42 = -0.9 , b43 = 1.2 ,

          b51 = -11.0/54.0 , b52 = 2.5 , b53 = -70.0/27.0 ,
          b54 = 35.0/27.0 ,

          b61 = 1631.0/55296.0 , b62 =   175.0/512.0 ,
          b63 =  575.0/13824.0 , b64 = 44275.0/110592.0 ,
          b65 =  253.0/4096.0 ,

          c1 = 37.0/378.0 , c3 = 250.0/621.0 , c4 = 125.0/594.0 ,
          c6 = 512.0/1771.0 ,
          dc5 = -277.0/14336.0 ;

    const double dc1 = c1 - 2825.0/27648.0 ,  
          dc3 = c3 - 18575.0/48384.0 ,
          dc4 = c4 - 13525.0/55296.0 , 
          dc6 = c6 - 0.25 ;

    // Initialise time to t0, needed when it is not updated by the integration.
    //       [ Note: Only for time dependent fields (usually electric) 
    //                 is it neccessary to integrate the time.] 
    //yOut[7] = yTemp[7]   = yIn[7]; 

    //  Saving yInput because yInput and yOut can be aliases for same array
    for(i=0;i<Nvar;i++) 
    {
        yIn[i]=yInput[i];
    }
    // RightHandSideInl(yIn, dydx) ;              // 1st Step

    for(i=0;i<Nvar;i++) 
    {
        yTemp[i] = yIn[i] + b21*Step*dydx[i] ;
    }
    this->RightHandSideInl(yTemp, ak2) ;              // 2nd Step

    for(i=0;i<Nvar;i++)
    {
        yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
    }
    this->RightHandSideInl(yTemp, ak3) ;              // 3rd Step

    for(i=0;i<Nvar;i++)
    {
        yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
    }
    this->RightHandSideInl(yTemp, ak4) ;              // 4th Step

    for(i=0;i<Nvar;i++)
    {
        yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                b54*ak4[i]) ;
    }
    this->RightHandSideInl(yTemp, ak5) ;              // 5th Step

    for(i=0;i<Nvar;i++)
    {
        yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                b64*ak4[i] + b65*ak5[i]) ;
    }
    this->RightHandSideInl(yTemp, ak6) ;              // 6th Step

    for(i=0;i<Nvar;i++)
    {
        // Accumulate increments with proper weights

        yOut[i] = yIn[i] + Step*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]) ;
    }
    for(i=0;i<Nvar;i++)
    {
        // Estimate error as difference between 4th and
        // 5th order methods

        yErr[i] = Step*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] +
                dc5*ak5[i] + dc6*ak6[i]) ;
    }
    for(i=0;i<Nvar;i++)
    {
        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = dydx[i];
    }
    fLastStepLength =Step;

    return ;
}

template <class T_Equation, int Nvar>
inline double
GUTCashKarpRKF45<T_Equation,Nvar>::
  DistChord()   const
{
    double distLine, distChord; 
    ThreeVector initialPoint, finalPoint, midPoint;

    // Store last initial and final points (they will be overwritten in self-Stepper call!)
    initialPoint = ThreeVector( fLastInitialVector[0], 
            fLastInitialVector[1], fLastInitialVector[2]); 
    finalPoint   = ThreeVector( fLastFinalVector[0],  
            fLastFinalVector[1],  fLastFinalVector[2]); 

    // Do half a step using StepNoErr

    fAuxStepper->GUTCashKarpRKF45::StepWithErrorEstimate( 
            fLastInitialVector, 
            fLastDyDx, 
            0.5 * fLastStepLength, 
            fMidVector,   
            fMidError );

    midPoint = ThreeVector( fMidVector[0], fMidVector[1], fMidVector[2]);       

    // Use stored values of Initial and Endpoint + new Midpoint to evaluate
    //  distance of Chord


    if (initialPoint != finalPoint) 
    {
        distLine  = GULineSection::Distline( midPoint, initialPoint, finalPoint );
        distChord = distLine;
    }
    else
    {
        distChord = (midPoint-initialPoint).Mag2();
    }
    return distChord;
}


#endif /*TCashKARP_RKF45 */
