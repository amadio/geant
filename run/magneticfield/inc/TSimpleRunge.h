#ifndef TSimpleRunge_HH
#define TSimpleRunge_HH

#include "TMagErrorStepper.h"
// #include "ThreeVector.h"

// #define  INTEGRATOR_CORRECTION   (1./((1<<2)-1))

template
<class T_Equation, unsigned int Nvar>
class TSimpleRunge : public TMagErrorStepper            
                     <TSimpleRunge<T_Equation, Nvar>, T_Equation, Nvar>
{
   public:  // with description
      static constexpr unsigned int OrderSimpleR= 2;
      static const     unsigned int Nmax_SR= 12;

      TSimpleRunge(T_Equation* EqRhs, unsigned int numStateVar = 0);
      TSimpleRunge(const TSimpleRunge& right);
      virtual ~TSimpleRunge(){ delete fEquation_Rhs; }

      virtual  VScalarIntegrationStepper* Clone() const;

      void SetEquationOfMotion(T_Equation* equation);

      inline  double IntegratorCorrection() { return  1./((1<<OrderSimpleR)-1); }
        
      inline __attribute__((always_inline)) 
      void RightHandSide(double y[], double charge, double dydx[]) 
      { fEquation_Rhs->T_Equation::RightHandSide(y, charge, dydx); }

      inline __attribute__((always_inline)) 
      void StepWithoutErrorEst( const double  yIn[],
                                      double charge,
                                const double  dydx[],
                                      double  h,
                                      double  yOut[]);

    private:
        //  Invariant(s) --- unchanged during simulation
        // Parameters
        unsigned int fNumberOfStateVariables;
        //  Owned objects - responsible for deleting!      
        T_Equation *fEquation_Rhs;

        //  State 
        double    yTemp[ Nvar>Nmax_SR ? Nvar : Nmax_SR ];
        double dydxTemp[ Nvar>Nmax_SR ? Nvar : Nmax_SR ];
           // scratch space    
};

//  Constructors

template <class T_Equation, unsigned int Nvar>
   TSimpleRunge<T_Equation,Nvar>::
   TSimpleRunge(T_Equation* EqRhs, unsigned int numStateVar)
   :
   TMagErrorStepper<TSimpleRunge<T_Equation, Nvar>, T_Equation, Nvar>
      (EqRhs, OrderSimpleR, (numStateVar > 0 ? numStateVar : Nvar) ),
   fNumberOfStateVariables(  numStateVar > 0 ? numStateVar : Nvar),
   fEquation_Rhs(EqRhs)
{
   //default GetNumberOfStateVariables() == Nmax_SR 
   assert (this->GetNumberOfStateVariables() <= Nmax_SR);
}

//  Copy constructor

template <class T_Equation, unsigned int Nvar>
   TSimpleRunge<T_Equation,Nvar>::
   TSimpleRunge(const TSimpleRunge& right)
   : TMagErrorStepper<TSimpleRunge<T_Equation, Nvar>, T_Equation, Nvar>( (T_Equation*) 0, 
                                                                         OrderSimpleR,
                                                                         right.fNumberOfStateVariables),
     fEquation_Rhs(new T_Equation(*(right.fEquation_Rhs)) )
{
   // Propagate it to the base class 
   TMagErrorStepper<TSimpleRunge<T_Equation, Nvar>, T_Equation, Nvar>
      ::SetEquationOfMotion(fEquation_Rhs);
   
   // SetEquationOfMotion(fEquation_Rhs);
}

template <class T_Equation, unsigned int Nvar>
   void TSimpleRunge<T_Equation,Nvar>::
     SetEquationOfMotion(T_Equation* equation)
{
   fEquation_Rhs= equation;
   TMagErrorStepper<TSimpleRunge<T_Equation, Nvar>, T_Equation, Nvar>
        ::SetEquationOfMotion(fEquation_Rhs);
   
   // TMagErrorStepper::SetEquationOfMotion(fEquation_Rhs);
}

template <class T_Equation, unsigned int Nvar>
VScalarIntegrationStepper* 
   TSimpleRunge<T_Equation,Nvar>::Clone() const
{
   return new TSimpleRunge<T_Equation,Nvar>( *this ); 
}

template <class T_Equation, unsigned int Nvar>
   inline __attribute__((always_inline)) 
   void TSimpleRunge<T_Equation,Nvar>::
        StepWithoutErrorEst( const double  yIn[],
                     double charge,
                     const double  dydx[],
                     double  h,
                     double  yOut[])
{
   // Initialise time to t0, needed when it is not updated by the integration.
   yTemp[7] = yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO
   
   for( unsigned int i = 0; i < Nvar; i++ )
   {
      yTemp[i] = yIn[i] + 0.5 * h*dydx[i] ;
   }
   this->RightHandSide(yTemp, charge, dydxTemp);
   
   for(unsigned int i = 0; i < Nvar; i++ )
   {
      yOut[i] = yIn[i] + h * ( dydxTemp[i] );
   }
}    
#endif /* TSimpleRunge_HH */
