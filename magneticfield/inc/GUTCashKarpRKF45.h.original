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

#include <iomanip>  // For  C++ style output (debug)
#include <iostream>

#include "GULineSection.h"
#include "GUVIntegrationStepper.h"
#include "base/Vector3D.h"

#define INLINERHS 1

#ifdef INLINERHS
#define REALLY_INLINE   inline __attribute__((always_inline))
#else
#define REALLY_INLINE   inline
#endif

template
<class T_Equation, unsigned int Nvar>
class GUTCashKarpRKF45 : public GUVIntegrationStepper
{
  public:
    static constexpr unsigned int sOrderMethod= 4;
    static constexpr unsigned int sNstore = (GUIntegrationNms::NumVarBase > Nvar) ? GUIntegrationNms::NumVarBase : Nvar;
                        // std::max( GUIntegrationNms::NumVarBase,  Nvar);
    // static const IntegratorCorrection = 1./((1<<4)-1);
    inline double IntegratorCorrection() { return 1./((1<<4)-1); }
    inline void   SetVerbose(bool v) { fVerbose= v; }
    
  public:
    inline
    GUTCashKarpRKF45( T_Equation *EqRhs,
                      unsigned int numStateVariables=0,
                      bool verbose=false);

    GUTCashKarpRKF45( const GUTCashKarpRKF45& );
    
    virtual ~GUTCashKarpRKF45();

    GUVIntegrationStepper* Clone() const override;

    REALLY_INLINE
       void StepWithErrorEstimate(const double* yInput,    // Consider __restrict__
                                  const double*  dydx,
                                        double   Step,
                                        double*  yOut,
                                        double*  yErr) override;

    double  DistChord()   const override;  

    REALLY_INLINE
    void RightHandSideInl(const double y[], double dydx[]) 
    {fEquation_Rhs->T_Equation::RightHandSide(y, dydx);}

    REALLY_INLINE
    void RightHandSideInl(const double y[], double dydx[], vecgeom::Vector3D<float> &Bfield) 
    { 
       // vecgeom::Vector3D<double> Position = { y[0], y[1], y[2] } ;
       PositionTmp.Set( y[0], y[1], y[2] );
       fEquation_Rhs->T_Equation::RightHandSide(y, PositionTmp, dydx, Bfield);
       // fEquation_Rhs->T_Equation::RightHandSide(y, dydx, Bfield);       
       // fEquation_Rhs->GetField()->T_Field::GetFieldValue(Point, Bfield);
       // fEquation_Rhs->TEvaluateRhsGivenB( y, Bfield, dydx ); 
    }
    
    void SetEquationOfMotion(T_Equation* equation);

    void PrintField( const char* label, const double y[6], const vecgeom::Vector3D<float> &Bfield ) const;
    void PrintDyDx( const char* label, const double dydx[Nvar], const double y[Nvar] ) const;
    void PrintDyDxLong( const char* label, const double dydx[Nvar], const double y[Nvar] ) const;     
    
    private:

    
      GUTCashKarpRKF45& operator=(const GUTCashKarpRKF45&) = delete;
        //private assignment operator.

    private:
        // 'Invariant' during integration - the pointers must not change
        // -----------
        T_Equation* fEquation_Rhs;
        bool        fOwnTheEquation;
        GUTCashKarpRKF45* fAuxStepper;

        // State -- intermediate values used during RK step
        // -----        
        double ak2[sNstore];
        double ak3[sNstore];
        double ak4[sNstore];
        double ak5[sNstore];
        double ak6[sNstore];
        double ak7[sNstore];
        double yTemp2[sNstore];
        double yTemp3[sNstore];
        double yTemp4[sNstore];
        double yTemp5[sNstore];
        double yTemp6[sNstore];
        double yIn[sNstore];
        vecgeom::Vector3D<float>  Bfield2, Bfield3, Bfield4, Bfield5, Bfield6;
        vecgeom::Vector3D<double>  PositionTmp;
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

        // Parameters - for debugging etc
        bool    fVerbose;
};

template <class T_Equation, unsigned int Nvar>
inline
GUTCashKarpRKF45<T_Equation,Nvar>::
   GUTCashKarpRKF45( T_Equation *EqRhs,
                     // unsigned int noIntegrationVariables,
                     unsigned int numStateVariables,
                     bool verbose)
   : GUVIntegrationStepper( EqRhs,    // dynamic_cast<GUVEquationOfMotion*>(EqRhs),
                            sOrderMethod,
                            Nvar,
                            ((numStateVariables>0) ? numStateVariables : sNstore) ),
     fEquation_Rhs(EqRhs),
     fOwnTheEquation(true),
     fAuxStepper(0),
     fLastStepLength(0.),
     fVerbose(verbose)
{
   assert( dynamic_cast<GUVEquationOfMotion*>(EqRhs) != 0 );
   assert( (numStateVariables == 0) || (numStateVariables >= Nvar) );
      
   fLastInitialVector = new double[sNstore] ;
   fLastFinalVector = new double[sNstore] ;
   fLastDyDx = new double[sNstore];
   
   fMidVector = new double[sNstore];
   fMidError =  new double[sNstore];
   if( verbose )
      std::cout << " GUTCashKarpRKF45 - constructed class. " << std::endl
                << " Nvar = " << Nvar << " Nstore= " << sNstore 
                << std::endl;
}

template <class T_Equation, unsigned int Nvar>
   void GUTCashKarpRKF45<T_Equation,Nvar>::
     SetEquationOfMotion(T_Equation* equation)
{
   fEquation_Rhs= equation;
   this->GUVIntegrationStepper::SetEquationOfMotion(fEquation_Rhs);
}

//  Copy - Constructor
// 
template <class T_Equation,unsigned int Nvar>
inline
GUTCashKarpRKF45<T_Equation,Nvar>::
   GUTCashKarpRKF45( const GUTCashKarpRKF45& right )
   : GUVIntegrationStepper( (GUVEquationOfMotion*) 0,
                            sOrderMethod,
                            Nvar,
                            right.GetNumberOfStateVariables() ),
     fEquation_Rhs( (T_Equation*) 0 ),
     fOwnTheEquation(true),
     fAuxStepper(0),   //  May overwrite below
     fLastStepLength(0.),
     fVerbose(right.fVerbose)
{
   SetEquationOfMotion( new T_Equation( *(right.fEquation_Rhs)) );
   fOwnTheEquation=true;
    // fEquation_Rhs= right.GetEquationOfMotion()->Clone());
   
   assert( dynamic_cast<GUVEquationOfMotion*>(fEquation_Rhs) != 0 );  
   assert( GetNumberOfStateVariables() >= Nvar);
      
   fLastInitialVector = new double[sNstore] ;
   fLastFinalVector = new double[sNstore] ;
   fLastDyDx = new double[sNstore];
   
   fMidVector = new double[sNstore];
   fMidError =  new double[sNstore];

   if( fVerbose )
      std::cout << " GUTCashKarpRKF45 - copy constructor: " << std::endl
                << " Nvar = " << Nvar << " Nstore= " << sNstore 
                << " Own-the-Equation = " << fOwnTheEquation << std::endl;

   if( right.fAuxStepper )
   {
      // Reuse the Equation of motion in the Auxiliary Stepper
      fAuxStepper = new GUTCashKarpRKF45(fEquation_Rhs, GetNumberOfStateVariables(), false);
   }
}

template <class T_Equation,unsigned int Nvar>
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
   if( fOwnTheEquation )
      delete fEquation_Rhs; // Expect to own the equation, except if auxiliary (then sharing the equation)
}

template <class T_Equation, unsigned int Nvar>
GUVIntegrationStepper* 
   GUTCashKarpRKF45<T_Equation,Nvar>::Clone() const
{
   // return new GUTCashKarpRKF45( *this );
   return new GUTCashKarpRKF45<T_Equation,Nvar>( *this );   
}


template <class T_Equation, unsigned int Nvar>
inline void
GUTCashKarpRKF45<T_Equation,Nvar>::
   StepWithErrorEstimate(const double*  yInput, // [],    
                         const double*  dydx, // [],
                         double Step,
                         double*  yOut, // [],
                         double*  yErr) // [])
{
    // const int nvar = 6 ;
    // const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
    std::cout << " Entered StepWithErrorEstimate of scalar " << std::endl;

    unsigned int i;

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

    std::cout<< " constants declared " <<std::endl;

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
    // PrintDyDx("dydx", dydx, yIn);
    std::cout<< " yin made " << std::endl;
#if 0
    double ak1[sNstore];
    vecgeom::Vector3D<float>  Bfield1;
    this->RightHandSideInl(yIn, ak1, Bfield1 );   // -- Get it again, for debugging
    // PrintField("yIn   ", yIn, Bfield1);
    // PrintDyDx("ak1-", ak1, yIn);
#endif
    
    std::cout<<" empty if else " << std::endl;
    for(i=0;i<Nvar;i++) 
    {
        yTemp2[i] = yIn[i] + b21*Step*dydx[i] ;
    }
    std::cout<<" just before rhs calculation " << std::endl;
    this->RightHandSideInl(yTemp2, ak2, Bfield2 ); // 2nd Step
    // PrintField("yTemp2", yTemp2, Bfield2);
    // PrintDyDx("ak2", ak2, yTemp2); 
    std::cout<<" 1 RHS calculating " << std::endl;
    for(i=0;i<Nvar;i++)
    {
        yTemp3[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
    }
    this->RightHandSideInl(yTemp3, ak3, Bfield3 ); // 3rd Step
    // PrintField("yTemp3", yTemp3, Bfield3);
    // PrintDyDx("ak3", ak3, yTemp3); 
    
    for(i=0;i<Nvar;i++)
    {
        yTemp4[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
    }
    this->RightHandSideInl(yTemp4, ak4, Bfield4 ); // 4th Step
    // PrintField("yTemp4", yTemp4, Bfield4);
    // PrintDyDx("ak4", ak4, yTemp4); 
    
    for(i=0;i<Nvar;i++)
    {
        yTemp5[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                b54*ak4[i]) ;
    }
    this->RightHandSideInl(yTemp5, ak5, Bfield5 );              // 5th Step
    // PrintField("yTemp5", yTemp5, Bfield5);
    // PrintDyDx("ak5", ak5, yTemp5); 

    for(i=0;i<Nvar;i++)
    {
        yTemp6[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                b64*ak4[i] + b65*ak5[i]) ;
    }
    this->RightHandSideInl(yTemp6, ak6, Bfield6 );              // 6th Step
    // PrintField("yTemp6", yTemp6, Bfield6);
    // PrintDyDx("ak6", ak6, yTemp6);

    for(i=0;i<Nvar;i++)
    {
        // Accumulate increments with proper weights

        yOut[i] = yIn[i] + Step * ( c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i] ) ;
    }
    for(i=0;i<Nvar;i++)
    {
        // Estimate error as difference between 4th and
        // 5th order methods

        yErr[i] = Step*( dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] +
                dc5*ak5[i] + dc6*ak6[i] ) ;
    }
    for(i=0;i<Nvar;i++)
    {
        // Store Input and Final values, for possible use in calculating chord
        fLastInitialVector[i] = yIn[i] ;
        fLastFinalVector[i]   = yOut[i];
        fLastDyDx[i]          = dydx[i];
    }
    fLastStepLength =Step;
    std::cout << " Exiting StepWithErrorEstimate of scalar " << std::endl;

    return ;
}

template <class T_Equation, unsigned int Nvar>
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



template <class T_Equation, unsigned int Nvar>
inline void
GUTCashKarpRKF45<T_Equation,Nvar>::
  PrintField( const char* label, const double y[Nvar], const vecgeom::Vector3D<float> &Bfield ) const
{
   std::cout << " PrintField/Stepper>  Field " << label << " "
             << "at x,y,z= ( " << y[0] << " , " << y[1]<< " , " << y[2] << " ) "
             << " is ( " << Bfield.x() << " , " << Bfield.y() << " , " << Bfield.z()
             <<  " ) kGauss - mag = " << Bfield.Mag() << std::endl;
}

template <class T_Equation, unsigned int Nvar>
inline void
GUTCashKarpRKF45<T_Equation,Nvar>::
   PrintDyDx( const char* label, const double dydx[Nvar], const double y[Nvar] ) const
{
   using std::cout;

   if( fVerbose > 0 )
   {    
      vecgeom::Vector3D<double>  dir(  dydx[0], dydx[1], dydx[2] );
      vecgeom::Vector3D<double>  dpds( dydx[3], dydx[4], dydx[5] );
      vecgeom::Vector3D<double>  p( y[3], y[4], y[5] );
      int oldPrec= cout.precision(3);   
      cout << " DyDx " << std::setw(4) << label << "> "
                << " xyz: " << std::setw(12) << dydx[0] << " , "
                << std::setw(12) << dydx[1] << " , "
                << std::setw(12) << dydx[2] << " ) "
                << " - mag = " << std::setw(12) << dir.Mag()
                << " dp/ds: " << std::setw(12) << dydx[3] << " , "
                << std::setw(12) << dydx[4] << " , "
                << std::setw(12) << dydx[5] << " "
                <<  " - mag = " << std::setw(5) << dpds.Mag()
                << " p-mag= " << p.Mag() << std::endl;
      cout.precision(oldPrec);
   }
}

template <class T_Equation, unsigned int Nvar>
inline void
GUTCashKarpRKF45<T_Equation,Nvar>::
   PrintDyDxLong( const char* label, const double dydx[Nvar], const double y[Nvar] ) const
{
   vecgeom::Vector3D<double>  dir(  dydx[0], dydx[1], dydx[2] );
   vecgeom::Vector3D<double>  dpds( dydx[3], dydx[4], dydx[5] );    
   std::cout << " PrintDyDx/Stepper>  dy/dx '" << std::setw(4) << label << "' "
             << " for x,y,z= ( " << dydx[0] << " , " << dydx[1]<< " , " << dydx[2] << " ) "
             << " - mag = " << dir.Mag() << std::endl
             << "                              "
             << " dp/ds(x,y,z) = ( " << dydx[3] << " , " << dydx[4]<< " , " << dydx[5] << " ) "
             <<  " ) - mag = " << dpds.Mag() << std::endl;
   vecgeom::Vector3D<double>  p( y[3], y[4], y[5] );
   double pMag= p.Mag();
   if( pMag == 0.0 ) pMag= 1.0;
   std::cout << "                                 "   
             << " 1/p dp/ds = " << dydx[3]/pMag << " , " << dydx[4]/pMag << " , " << dydx[5]/pMag << " ) "
             << std::endl;
   std::cout << "                                 "      
             << "         p = " << y[3] << " , " << y[4] << " , " << y[5] << " ) " << std::endl;
}
#endif /*TCashKARP_RKF45 */
