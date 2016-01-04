#ifndef TemplateTSimpleRunge_HH
#define TemplateTSimpleRunge_HH

#include "TemplateTMagErrorStepper.h"

// #define  INTEGRATOR_CORRECTION   (1./((1<<2)-1))

template
<class Backend, class T_Equation, unsigned int Nvar>
class TemplateTSimpleRunge : public TemplateTMagErrorStepper            
                     <Backend, TemplateTSimpleRunge<Backend, T_Equation, Nvar>, T_Equation, Nvar>
{
   public:  // with description

      typedef typename Backend::precision_v  Double_v;

      static constexpr unsigned int OrderSimpleR= 2;
      static const     unsigned int Nmax_SR= 12;

      TemplateTSimpleRunge(T_Equation* EqRhs, unsigned int numStateVar = 0);
      TemplateTSimpleRunge(const TemplateTSimpleRunge& right);
      virtual ~TemplateTSimpleRunge(){ delete fEquation_Rhs; }

      virtual  TemplateGUVIntegrationStepper<Backend>* Clone() const;

      void SetEquationOfMotion(T_Equation* equation);

      inline  double IntegratorCorrection() { return  1./((1<<OrderSimpleR)-1); }
        
      inline __attribute__((always_inline)) 
      void RightHandSide(Double_v y[],/* Double_v charge,*/ Double_v dydx[]) 
      { fEquation_Rhs->T_Equation::RightHandSide(y,/* charge,*/ dydx); }

      inline __attribute__((always_inline)) 
      void StepWithoutErrorEst( const Double_v  yIn[],
                                const Double_v  dydx[],
                                      Double_v  h,
                                      Double_v  yOut[]);

    private:
        //  Invariant(s) --- unchanged during simulation
        // Parameters
        unsigned int fNumberOfStateVariables;
        //  Owned objects - responsible for deleting!      
        T_Equation *fEquation_Rhs;

        //  State 
        Double_v    yTemp[ Nvar>Nmax_SR ? Nvar : Nmax_SR ];
        Double_v dydxTemp[ Nvar>Nmax_SR ? Nvar : Nmax_SR ];
           // scratch space    
};

//  Constructors

template <class Backend, class T_Equation, unsigned int Nvar>
TemplateTSimpleRunge<Backend,T_Equation,Nvar>
  ::TemplateTSimpleRunge(T_Equation* EqRhs, unsigned int numStateVar)
   : TemplateTMagErrorStepper<Backend, TemplateTSimpleRunge<Backend, T_Equation, Nvar>, T_Equation, Nvar>
      (EqRhs, OrderSimpleR, (numStateVar > 0 ? numStateVar : Nvar) ),
     fNumberOfStateVariables(  numStateVar > 0 ? numStateVar : Nvar),
     fEquation_Rhs(EqRhs)
{
   //default GetNumberOfStateVariables() == Nmax_SR 
   assert (this->GetNumberOfStateVariables() <= Nmax_SR);
}

//  Copy constructor

template <class Backend, class T_Equation, unsigned int Nvar>
TemplateTSimpleRunge<Backend,T_Equation,Nvar>
  ::TemplateTSimpleRunge(const TemplateTSimpleRunge& right)
   : TemplateTMagErrorStepper<Backend, TemplateTSimpleRunge<Backend, T_Equation, Nvar>, T_Equation, Nvar>
      ( (T_Equation*) 0, 
        OrderSimpleR,
        right.fNumberOfStateVariables),
     fEquation_Rhs(new T_Equation(*(right.fEquation_Rhs)) )
{
   // Propagate it to the base class 
   TemplateTMagErrorStepper<Backend, TemplateTSimpleRunge<Backend, T_Equation, Nvar>, T_Equation, Nvar>
      ::SetEquationOfMotion(fEquation_Rhs);   
}

template <class Backend, class T_Equation, unsigned int Nvar>
void TemplateTSimpleRunge<Backend,T_Equation,Nvar>
  ::SetEquationOfMotion(T_Equation* equation)
{
   fEquation_Rhs= equation;
   TemplateTMagErrorStepper<Backend, TemplateTSimpleRunge<Backend, T_Equation, Nvar>, T_Equation, Nvar>
        ::SetEquationOfMotion(fEquation_Rhs);
}

template <class Backend, class T_Equation, unsigned int Nvar>
TemplateGUVIntegrationStepper<Backend>* 
   TemplateTSimpleRunge<Backend,T_Equation,Nvar>::Clone() const
{
   return new TemplateTSimpleRunge<Backend,T_Equation,Nvar>( *this ); 
}

template <class Backend, class T_Equation, unsigned int Nvar>
inline __attribute__((always_inline)) 
void TemplateTSimpleRunge<Backend,T_Equation,Nvar>
  ::StepWithoutErrorEst( const typename Backend::precision_v  yIn[],
                         const typename Backend::precision_v  dydx[],
                               typename Backend::precision_v  h,
                               typename Backend::precision_v  yOut[])
{
   // Initialise time to t0, needed when it is not updated by the integration.
   yTemp[7] = yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO
   
   for( unsigned int i = 0; i < Nvar; i++ )
   {
      yTemp[i] = yIn[i] + 0.5 * h*dydx[i] ;
   }
   this->RightHandSide(yTemp,dydxTemp);
   
   for(unsigned int i = 0; i < Nvar; i++ )
   {
      yOut[i] = yIn[i] + h * ( dydxTemp[i] );
   }
}    
#endif /* TemplateTSimpleRunge_HH */
