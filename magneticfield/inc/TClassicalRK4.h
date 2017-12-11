#ifndef TCLASSICALRK4_H
#define TCLASSICALRK4_H

#include "Geant/Config.h"
#include "TMagErrorStepper.h"

// #include <algorithm> // for std::max

// #define  INTEGRATOR_CORRECTION   (1./((1<<2)-1))

// template <class T> inline constexpr const T& MaxConst (const T& a, const T& b) { return (a<b)?b:a;  } 

template
<class T_Equation, unsigned int Nvar>
class TClassicalRK4 : public  TMagErrorStepper
                      <TClassicalRK4<T_Equation, Nvar>, T_Equation, Nvar>
{
  public:  // with description
    static constexpr unsigned int OrderRK4= 4;
    static constexpr unsigned int NumVarStore = (GUIntegrationNms::NumVarBase > Nvar) ? GUIntegrationNms::NumVarBase : Nvar;
                        // MaxConst( GUIntegrationNms::NumVarBase,  Nvar);
                        // std::max( GUIntegrationNms::NumVarBase,  Nvar);

    TClassicalRK4(T_Equation *EqRhs) // , int numberOfVariables = 8)
       : TMagErrorStepper<TClassicalRK4<T_Equation, Nvar>, T_Equation, Nvar>(EqRhs, OrderRK4, Nvar)
       // fEquation_Rhs(EqRhs)
    {
    }

    TClassicalRK4(const TClassicalRK4& right);

    virtual  VScalarIntegrationStepper* Clone() const override final;
    
    // void SetOurEquationOfMotion(T_Equation* equation);
       
    virtual ~TClassicalRK4(){ } // delete fEquation_Rhs;}

    // static const IntegratorCorrection = 1./((1<<4)-1); 
    inline double IntegratorCorrection() { return 1./((1<<OrderRK4)-1); }
    
    // A stepper that does not know about errors.
    // It is used by the MagErrorStepper stepper.
    inline
    void  StepWithoutErrorEst( const double  yIn[],
                               double charge,
                               const double  dydx[],
                               double  h,
                               double  yOut[]);  // override final;  => Not virtual method, must exist though!
 
  public:
    // __attribute__((always_inline)) 
    //  int IntegratorOrder() const { return OrderRK4; }

  private:
        TClassicalRK4& operator=(const TClassicalRK4&) = delete;
        // Private assignment operator.

    private:
        // Invariants
        static constexpr unsigned int Nvarstor= 8 * ((Nvar-1)/8+1); 
        
        // Owned Object
        //  T_Equation *fEquation_Rhs;

        // STATE
        
        // scratch space
        double dydxm[Nvarstor]; 
        double dydxt[Nvarstor]; 
        double yt[Nvarstor];
};

/*
template <class T_Equation, unsigned int Nvar>
   void TClassicalRK4<T_Equation,Nvar>::
     SetOurEquationOfMotion(T_Equation* equation)
{
   // fEquation_Rhs= equation;
   TMagErrorStepper<TClassicalRK4<T_Equation, Nvar>, T_Equation, Nvar>
       ::SetEquationOfMotion(equation); // fEquation_Rhs);
   
   // TMagErrorStepper::SetEquationOfMotion(fEquation_Rhs);
}
*/

template <class T_Equation, unsigned int Nvar>
  TClassicalRK4<T_Equation,Nvar>::
  TClassicalRK4(const TClassicalRK4& right)
   :
      TMagErrorStepper<TClassicalRK4<T_Equation, Nvar>, T_Equation, Nvar>( // (T_Equation*) 0,
           new T_Equation(*(right.fEquation_Rhs)),
           OrderRK4,
           right.GetNumberOfStateVariables() )
   // TMagErrorStepper<TClassicalRK4<T_Equation, Nvar>, T_Equation, Nvar>( right ), 

   //  right.fEquation_Rhs->Clone())  
   //  right.fEquation_Rhs->CloneOrSafeSelf())  // Alt   
{
   // TMagErrorStepper<TClassicalRK4<T_Equation, Nvar>, T_Equation, Nvar>
   //   ::SetEquationOfMotion(fEquation_Rhs);
   // assert(fEquation_Rhs);
   // SetEquationOfMotion(fEquation_Rhs);    // Propagates it also to base class 

   // TMagErrorStepper::SetEquationOfMotion(fEquation_Rhs);
}  

template <class T_Equation, unsigned int Nvar>
VScalarIntegrationStepper* 
   TClassicalRK4<T_Equation,Nvar>::Clone() const
{
   // return new TClassicalRK4<T_Equation,Nvar>( *this );
   auto clone= new TClassicalRK4<T_Equation,Nvar>( *this );
   // clone->Check();
   assert ( clone->fEquation_Rhs != 0 );
   return clone;
}

static constexpr double inv6=1./6;

#define INLINEDUMBSTEPPER 1

template <class T_Equation, unsigned int Nvar>
#ifdef INLINEDUMBSTEPPER
   GEANT_FORCE_INLINE
#else
// __attribute__((noinline))
#endif 
   void TClassicalRK4<T_Equation,Nvar>::StepWithoutErrorEst( const double  yIn[],
                                                  double charge,
                                                  const double  dydx[],
                                                  double  h,
                                                  double  yOut[])
   // Given values for the variables y[0,..,n-1] and their derivatives
   // dydx[0,...,n-1] known at x, use the classical 4th Runge-Kutta
   // method to advance the solution over an interval h and return the
   // incremented variables as yout[0,...,n-1], which not be a distinct
   // array from y. The user supplies the routine RightHandSide(x,y,dydx),
   // which returns derivatives dydx at x. The source is routine rk4 from
   // NRC p. 712-713 .
{
   unsigned int i;
   constexpr double one_sixth= 1.0 / 6.0; 
   double  hh = h*0.5 , h6 = h * one_sixth; // (1.0/6.0);

   // Initialise time to t0, needed when it is not updated by the integration.
   //        [ Note: Only for time dependent fields (usually electric) 
   //                  is it neccessary to integrate the time.] 
   // yt[7]   = yIn[7]; 
   //yOut[7] = yIn[7];
   
   for(i=0;i<Nvar;i++)
   {
      yt[i] = yIn[i] + hh*dydx[i] ;             // 1st Step K1=h*dydx
   }
   this->RightHandSide(yt, charge, dydxt) ;              // 2nd Step K2=h*dydxt
   
   for(i=0;i<Nvar;i++)
   { 
      yt[i] = yIn[i] + hh*dydxt[i] ;
   }
   this->RightHandSide(yt, charge, dydxm) ;              // 3rd Step K3=h*dydxm

   for(i=0;i<Nvar;i++)
   {
      yt[i]   = yIn[i] + h*dydxm[i] ;
      dydxm[i] += dydxt[i] ;                    // now dydxm=(K2+K3)/h
   }
   this->RightHandSide(yt, charge, dydxt) ;              // 4th Step K4=h*dydxt
   
   for(i=0;i<Nvar;i++)    // Final RK4 output
   {
      yOut[i] = yIn[i]+h6*(dydx[i]+dydxt[i]+2.0*dydxm[i]); //+K1/6+K4/6+(K2+K3)/3
   }
   //            if ( Nvar == 12 )  { this->NormalisePolarizationVector ( yOut ); }
   
}  // end of DumbStepper ....................................................



#endif
