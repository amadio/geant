//
// Embedded explicit Runge-Kutta Stepper using Cash Karp's RK tableau
//
// Adapted to remove 'backend' template interface: J. Apostolakis, Oct 2017
// Templated version:  Ananya, Feb/March 2016
//     ( commit 95e1316bcc156a04c876d6ea0fc9e60a15eeac4f )
// 
// Adapted from 'GUTCashKarpRKF45' by John Apostolakis, Nov 2015
//
// Adapted from 'GUTCashKarpRKF45' by Qieshen Xie, GSoC 2014
//         (derived from G4CashKarpRKF45)
//
//     
#ifndef VECTORCASHKARPRKF45_H
#define VECTORCASHKARPRKF45_H

#include "GUVectorLineSection.h"
#include "GUVVectorIntegrationStepper.h"

// #include "TMagErrorStepper.h" //for sake of GUIntegrationNms::NumVars
// #include "TemplateTMagErrorStepper.h"

// #include "AlignedBase.h"  // ==> Ensures alignment of storage for Vc objects

template
<class T_Equation, unsigned int Nvar>
   class VectorCashKarpRKF45 : public GUVVectorIntegrationStepper // <Backend>, public AlignedBase
{
  public:
   
    template <typename T>
       using Vector3D = vecgeom::Vector3D<T>;
  
    using Double_v        = Geant::Double_v;
    using ThreeVectorSimd = Vector3D<Double_v>;
  
    static constexpr unsigned int sOrderMethod= 4;
    static constexpr unsigned int sNstore = 6;   // How many variables the full state entails
               // (GUIntegrationNms::NumVarBase > Nvar) ? GUIntegrationNms::NumVarBase : Nvar;
    // std::max( GUIntegrationNms::NumVarBase,  Nvar);
    // static const double IntegratorCorrection = 1./((1<<4)-1);
    inline double IntegratorCorrection() { return 1./((1<<sOrderMethod)-1); }

  public:
    inline
    VectorCashKarpRKF45( T_Equation *EqRhs,
                         unsigned int numStateVariables=0);

    VectorCashKarpRKF45( const VectorCashKarpRKF45& );
    
    virtual ~VectorCashKarpRKF45();

    GUVVectorIntegrationStepper* Clone() const override;

    GEANT_FORCE_INLINE
    void StepWithErrorEstimate(const Double_v  yInput[],  // Consider __restrict__
                               const Double_v  dydx[],
                               const Double_v& charge,
                               const Double_v& hStep,
                                     Double_v  yOut[],
                                     Double_v  yErr[]) override;

    Double_v DistChord() const override ;

    GEANT_FORCE_INLINE
    void RightHandSideInl(Double_v y[], const Double_v& charge, Double_v dydx[]) 
    {fEquation_Rhs->T_Equation::RightHandSide(y, charge, dydx);}

    void SetEquationOfMotion(T_Equation* equation);
    
    private:
      VectorCashKarpRKF45& operator=(const VectorCashKarpRKF45&) = delete;
        //private assignment operator.

    private: 
         // State -- intermediate values used during RK step
         // -----        
         Double_v ak2[sNstore];
         Double_v ak3[sNstore];
         Double_v ak4[sNstore];
         Double_v ak5[sNstore];
         Double_v ak6[sNstore];
         Double_v ak7[sNstore];
         Double_v yTemp2[sNstore];  // Separate temporaries per step - to aid compiler
         Double_v yTemp3[sNstore];
         Double_v yTemp4[sNstore];
         Double_v yTemp5[sNstore];
         Double_v yTemp6[sNstore];                
         
         Double_v yIn[sNstore];
         // scratch space
         
         // State -- values used for subsequent call to DistChord
         // -----
         
         Double_v  fLastStepLength;
         
         Double_v  fLastInitialVector[sNstore];
         Double_v  fLastFinalVector[sNstore];
         Double_v  fLastDyDx[sNstore];
         Double_v  fMidVector[sNstore];
         Double_v  fMidError[sNstore];
         // for DistChord calculations
      
    private:
        // 'Invariant' during integration - the pointers must not change
        // -----------
        T_Equation* fEquation_Rhs;
        bool        fOwnTheEquation;  //  --> indicates ownership of Equation object

        bool  fDebug= false;
};

template <class T_Equation, unsigned int Nvar>
inline
VectorCashKarpRKF45<T_Equation,Nvar>::
   VectorCashKarpRKF45( T_Equation *EqRhs,
                        unsigned int numStateVariables )
   : GUVVectorIntegrationStepper( nullptr, // EqRhs,   ==>>  Does not inherit !!
                                  sOrderMethod,
                                  Nvar, //8, //Ananya
                                  ((numStateVariables>0) ? numStateVariables : sNstore) ),
     fEquation_Rhs(EqRhs),
     // fLastStepLength(0.),
     fOwnTheEquation(false)
{
   if( fDebug ) {
      std::cout<<"\n----Entered constructor of VectorCashKarpRKF45 "<<std::endl;
      std::cout<<"----In VectorCashKarpRKF45 constructor, Nvar is: "<<Nvar<<std::endl;
   }
   // assert( dynamic_cast<TemplateGUVEquationOfMotion<Backend>*>(EqRhs) != 0 );

   fLastStepLength= Double_v(0.);
   
   assert( (numStateVariables == 0) || (numStateVariables >= Nvar) );
  
   std::cout<<"----end of constructor of VectorCashKarpRKF45"<<std::endl;
}

template <class T_Equation, unsigned int Nvar>
void VectorCashKarpRKF45<T_Equation,Nvar>::
  SetEquationOfMotion(T_Equation* equation)
{
   fEquation_Rhs= equation;
   this->GUVVectorIntegrationStepper::SetABCEquationOfMotion(nullptr); // fEquation_Rhs);
}

//  Copy - Constructor
// 
template <class T_Equation,unsigned int Nvar>
inline
VectorCashKarpRKF45<T_Equation,Nvar>::
   VectorCashKarpRKF45( const VectorCashKarpRKF45& right )
   : GUVVectorIntegrationStepper( (GUVVectorEquationOfMotion*) nullptr,
                                              sOrderMethod,
                                              Nvar,
                                              right.GetNumberOfStateVariables() ),
     fEquation_Rhs( (T_Equation*) nullptr ),
     fOwnTheEquation(false)
{
   if( fDebug ) {
    std::cout<<"----Entered constructor of VectorCashKarpRKF45 "<<std::endl;
   }
   SetEquationOfMotion( new T_Equation( *(right.fEquation_Rhs)) );
    // fEquation_Rhs= right.GetEquationOfMotion()->Clone());
   
   // assert( dynamic_cast<GUVVectorEquationOfMotion*>(fEquation_Rhs) != 0 );   // No longer Deriving
   assert( this->GetNumberOfStateVariables() >= Nvar);

   fLastStepLength= Double_v(0.);   
   
   if( fDebug )
      std::cout << " VectorCashKarpRKF45 - copy constructor: " << std::endl
                << " Nvar = " << Nvar << " Nstore= " << sNstore 
                << " Own-the-Equation = " << fOwnTheEquation << std::endl;
}



template <class T_Equation,unsigned int Nvar>
GEANT_FORCE_INLINE
VectorCashKarpRKF45<T_Equation,Nvar>::~VectorCashKarpRKF45()
{
   std::cout<<"----- Vector CashKarp destructor"<<std::endl;
   if( fOwnTheEquation )
      delete fEquation_Rhs; // Expect to own the equation, except if auxiliary (then sharing the equation)

  std::cout<<"----- VectorCashKarp destructor (ended)"<<std::endl;
}

template <class T_Equation, unsigned int Nvar>
GUVVectorIntegrationStepper* 
VectorCashKarpRKF45<T_Equation,Nvar>::Clone() const
{
   // return new VectorCashKarpRKF45( *this );
   return new VectorCashKarpRKF45<T_Equation,Nvar>( *this );   
}


template <class T_Equation, unsigned int Nvar>
inline void
VectorCashKarpRKF45<T_Equation,Nvar>::

  StepWithErrorEstimate(const Double_v  yInput[],    
                        const Double_v  dydx[],
                        const Double_v& charge,                        
                        const Double_v& Step,
                              Double_v  yOut[],
                              Double_v  yErr[])
{
    // const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
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

    // Initialise time to t0, needed when it is not updated by the integration.
    //       [ Note: Only for time dependent fields (usually electric) 
    //                 is it neccessary to integrate the time.] 
    //yOut[7] = yTemp[7]   = yIn[7]; 

    //  Saving yInput because yInput and yOut can be aliases for same array
    for(i=0;i<Nvar;i++) 
    {
        yIn[i]=yInput[i];
    }
    // RightHandSideInl(yIn, charge,  dydx) ;              // 1st Step

    for(i=0;i<Nvar;i++) 
    {
        yTemp2[i] = yIn[i] + b21*Step*dydx[i] ;
    }
    this->RightHandSideInl(yTemp2, charge,  ak2) ;              // 2nd Step

    for(i=0;i<Nvar;i++)
    {
        yTemp3[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
    }
    this->RightHandSideInl(yTemp3, charge,  ak3) ;              // 3rd Step

    for(i=0;i<Nvar;i++)
    {
        yTemp4[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
    }
    this->RightHandSideInl(yTemp4, charge,  ak4) ;              // 4th Step

    for(i=0;i<Nvar;i++)
    {
        yTemp5[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                b54*ak4[i]) ;
    }
    this->RightHandSideInl(yTemp5, charge,  ak5) ;              // 5th Step

    for(i=0;i<Nvar;i++)
    {
       yTemp6[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                                  b64*ak4[i] + b65*ak5[i]) ;
    }
    this->RightHandSideInl(yTemp6, charge,  ak6) ;              // 6th Step

    for(i=0;i<Nvar;i++)
    {
        // Accumulate increments with correct weights
        yOut[i] = yIn[i] + Step*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]) ;
    }
    for(i=0;i<Nvar;i++)
    {
        // Estimate error as difference between 4th and
        // 5th order methods

        yErr[i] = Step*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] +
                dc5*ak5[i] + dc6*ak6[i]) ;
        // std::cout<< "----In Stepper, yerrr is: "<<yErr[i]<<std::endl;
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

template <class T_Equation, unsigned int Nvar>
inline 
Geant::Double_v
VectorCashKarpRKF45<T_Equation,Nvar>::
  DistChord()  const
{
#if 1
    Double_v        distChord= 0.0;
#else   
    Double_v        distLine, distChord;
    ThreeVectorSimd initialPoint, finalPoint, midPoint;

    // Store last initial and final points (they will be overwritten in self-Stepper call!)
    initialPoint = ThreeVectorSimd( fLastInitialVector[0], 
                                    fLastInitialVector[1], fLastInitialVector[2]); 
    finalPoint   = ThreeVectorSimd( fLastFinalVector[0],  
                                    fLastFinalVector[1],  fLastFinalVector[2]); 

    // Do half a step using StepNoErr
    fAuxStepper->StepWithErrorEstimate( fLastInitialVector, 
                                        fLastDyDx, 
                                        0.5 * fLastStepLength, 
                                        fMidVector,   
                                        fMidError );

    midPoint = ThreeVectorSimd( fMidVector[0], fMidVector[1], fMidVector[2]);       

    // Use stored values of Initial and Endpoint + new Midpoint to evaluate
    //  distance of Chord

    distChord  = GUVectorLineSection::Distline( midPoint, initialPoint, finalPoint );
#endif    
    return distChord;
}

#endif /*GUV Vector CashKARP_RKF45 */
