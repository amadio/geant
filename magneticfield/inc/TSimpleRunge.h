#ifndef TSimpleRunge_HH
#define TSimpleRunge_HH

#include "TMagErrorStepper.h"
#include "ThreeVector.h"

#define  Nmax_SR   12
// #define  INTEGRATOR_CORRECTION   (1./((1<<2)-1))

template
<class T_Equation, int N>
class TSimpleRunge : public TMagErrorStepper            
                     <TSimpleRunge<T_Equation, N>, T_Equation, N>
{

    public:  // with description

        // static  const double IntegratorCorrection = 1./((1<<2)-1);
        inline  double IntegratorCorrection() { return  1./((1<<2)-1); } 
        // const   int    Nmax=12; 
        TSimpleRunge(T_Equation* EqRhs, 
                     int numberOfVariables = 6)
            :
              TMagErrorStepper
                 <TSimpleRunge<T_Equation, N>, T_Equation, N>
                  (EqRhs, numberOfVariables),
              fNumberOfVariables(numberOfVariables),
              fEquation_Rhs(EqRhs)
        {
            //default GetNumberOfStateVariables() == Nmax_SR 
            assert (this->GetNumberOfStateVariables() <= Nmax_SR);
        }


        ~TSimpleRunge(){;}


        __attribute__((always_inline)) 
        void 
        RightHandSide(double y[], double dydx[]) 
        { fEquation_Rhs->T_Equation::RightHandSide(y, dydx); }


        __attribute__((always_inline)) 
        void
        DumbStepper( const double  yIn[],
                     const double  dydx[],
                     double  h,
                     double  yOut[])
        {
            // Initialise time to t0, needed when it is not updated by the integration.
            yTemp[7] = yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO

            int i;

            for( i = 0; i < N; i++ ) 
            {
                yTemp[i] = yIn[i] + 0.5 * h*dydx[i] ;
            }

            this->RightHandSide(yTemp,dydxTemp);

            for( i = 0; i < N; i++ ) 
            {
                yOut[i] = yIn[i] + h * ( dydxTemp[i] );
            }
        } 

    public:  // without description

        __attribute__((always_inline)) 
        int  
        IntegratorOrder() const { return 2; }

    private:

        int fNumberOfVariables ;
        double dydxTemp[N>Nmax_SR?N:Nmax_SR];
        double    yTemp[N>Nmax_SR?N:Nmax_SR];

        T_Equation *fEquation_Rhs;
        // scratch space    
};

#endif /* TSimpleRunge_HH */
