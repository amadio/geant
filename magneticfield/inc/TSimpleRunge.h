#ifndef TSimpleRunge_HH
#define TSimpleRunge_HH

#include "TMagErrorStepper.h"
#include "ThreeVector.h"


// #define  INTEGRATOR_CORRECTION   (1./((1<<2)-1))

template
<class T_Equation, int Nvar>
class TSimpleRunge : public TMagErrorStepper            
                     <TSimpleRunge<T_Equation, Nvar>, T_Equation, Nvar>
{

    public:  // with description
        static constexpr unsigned int OrderSimpleR= 2;
        static const     unsigned int Nmax_SR= 12;
        // static constexpr double IntegratorCorrectionFactor = 1./((1<<OrderSimpleR)-1);
        inline  double IntegratorCorrection() { return  1./((1<<OrderSimpleR)-1); }
        // const   int    Nmax=12;
        
        TSimpleRunge(T_Equation* EqRhs, 
                     int numStateVar = -1)
            :
              TMagErrorStepper
                 <TSimpleRunge<T_Equation, Nvar>, T_Equation, Nvar>
                 (EqRhs, OrderSimpleR, Nvar, numStateVar > 0 ? numStateVar : Nvar ),
              fNumberOfStateVariables(numStateVar),
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
        StepWithoutErrorEst( const double  yIn[],
                     const double  dydx[],
                     double  h,
                     double  yOut[])
        {
            // Initialise time to t0, needed when it is not updated by the integration.
            yTemp[7] = yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO

            int i;

            for( i = 0; i < Nvar; i++ )
            {
                yTemp[i] = yIn[i] + 0.5 * h*dydx[i] ;
            }

            this->RightHandSide(yTemp,dydxTemp);

            for( i = 0; i < Nvar; i++ )
            {
                yOut[i] = yIn[i] + h * ( dydxTemp[i] );
            }
        } 

    public:  // without description

        __attribute__((always_inline)) 
        int  
        IntegratorOrder() const { return OrderSimpleR; }

    private:

        int fNumberOfStateVariables ;
        double dydxTemp[ Nvar>Nmax_SR ? Nvar : Nmax_SR ];
        double    yTemp[ Nvar>Nmax_SR ? Nvar : Nmax_SR ];

        T_Equation *fEquation_Rhs;
        // scratch space    
};

#endif /* TSimpleRunge_HH */
