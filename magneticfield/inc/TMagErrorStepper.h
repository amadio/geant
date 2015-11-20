#ifndef TMAGERRORSTEPPER_HH
#define TMAGERRORSTEPPER_HH

#include <algorithm> // for std::max

// #include "G4Types.hh"
#include "GUVIntegrationStepper.h"
#include "ThreeVector.h"
#include "GULineSection.h"

namespace GUIntegrationNms
{
   constexpr unsigned int NumVarBase  = 8;  //
}

template
<class T_Stepper, class T_Equation, unsigned int Nvar>
class TMagErrorStepper : public GUVIntegrationStepper
{
    public:  // with description
        static constexpr unsigned int NumVarStore = (Nvar > GUIntegrationNms::NumVarBase) ?
                                                     Nvar : GUIntegrationNms::NumVarBase ;
           // std::max( GUIntegrationNms::NumVarBase,  Nvar);

        TMagErrorStepper( T_Equation *EqRhs,
                          unsigned int integrationOrder,   // Make it a template Parameter ??
                          unsigned int numStateVariables); // = -1)  // No default -- must ensure order is set

        TMagErrorStepper( const TMagErrorStepper& right );
   
        virtual ~TMagErrorStepper() {;}

        inline void RightHandSide(double y[], double dydx[]) 
              { fEquation_Rhs->T_Equation::RightHandSide(y, dydx); }

        inline void StepWithErrorEstimate( const double yInput[],
                                           const double dydx[],
                                           double hstep,
                                           double yOutput[],
                                           double yError []      );
            // The stepper for the Runge Kutta integration. The stepsize 
            // is fixed, with the Step size given by h.
            // Integrates ODE starting values y[0 to 6].
            // Outputs yout[] and its estimated error yerr[].
        

        double DistChord() const; 

    private:
        TMagErrorStepper& operator=(const TMagErrorStepper&) = delete;
        // Private assignment operator.

    private:

        // STATE
        ThreeVector fInitialPoint, fMidPoint, fFinalPoint;
        // Data stored in order to find the chord

        // Dependent Objects, owned --- part of the STATE 
        double yInitial[NumVarStore];   // [Nvar<8?8:Nvar];
        double yMiddle[NumVarStore];
        double dydxMid[NumVarStore];
        double yOneStep[NumVarStore];
        // The following arrays are used only for temporary storage
        // they are allocated at the class level only for efficiency -
        // so that calls to new and delete are not made in Stepper().

        T_Equation *fEquation_Rhs;
};

template<class T_Stepper, class T_Equation, unsigned int Nvar>
   TMagErrorStepper<T_Stepper, T_Equation, Nvar>::
   TMagErrorStepper( T_Equation *EqRhs,
                     unsigned int integrationOrder,   // Make it a template Parameter ??
                     unsigned int numStateVariables) // = 0)  // No default -- must ensure order is set
   : GUVIntegrationStepper( EqRhs,
                            integrationOrder,
                            Nvar,                // Here we must pass it to base class !
                            numStateVariables ), // ((numStateVariables>0) ? numStateVariables : NumVarStore) ),
   fEquation_Rhs(EqRhs)
{
   assert( numStateVariables >= Nvar ); 
}

template<class T_Stepper, class T_Equation, unsigned int Nvar>
   TMagErrorStepper<T_Stepper, T_Equation, Nvar>::
   TMagErrorStepper( const TMagErrorStepper& right )
    :
       GUVIntegrationStepper( (T_Equation *) 0, 
                              right.IntegrationOrder(),
                              right.GetNumberOfVariables(),  // must be == Nvar
                              right.GetNumberOfStateVariables() ), 
       fEquation_Rhs(right.GetEquationOfMotion()->Clone())
{
   SetEquationOfMotion(fEquation_Rhs); 

   // unsigned nvar = std::max(this->GetNumberOfVariables(), 8);
   assert( this->GetNumberOfVariables() == Nvar ); 
}

// inline
template<class T_Stepper, class T_Equation, unsigned int Nvar>
void
   TMagErrorStepper<T_Stepper, T_Equation, Nvar>::
StepWithErrorEstimate( const double yInput[],
                const double dydx[],
                double hstep,
                double yOutput[],
                double yError []      )
            // The stepper for the Runge Kutta integration. The stepsize 
            // is fixed, with the Step size given by h.
            // Integrates ODE starting values y[0 to 6].
            // Outputs yout[] and its estimated error yerr[].
{  
   const unsigned maxvar= GetNumberOfStateVariables();

   // correction for Richardson Extrapolation.
   //double  correction = 1. / ( (1 << 
   //          static_cast<T_Stepper*>(this)->T_Stepper::IntegratorOrder()) -1 );
   //  Saving yInput because yInput and yOutput can be aliases for same array
   
   for(unsigned int i=0;i<NumVarStore;i++){
      yInitial[i]= yInput[i];
      yOutput[i] = yInput[i];
      yError[i]  = 0.0;         
   }
   
   // Copy the remaining state - part which is not integrated
   for(unsigned int i=Nvar+1;i<NumVarStore;i++){
      yMiddle[i]=yInput[i];   
      yOneStep[i] = yInput[i]; // As it contributes to final value of yOutput ?
   }

   // const unsigned maxvar= GetNumberOfStateVariables();
   // for(i=Nvar;i<maxvar;i++) yOutput[i]=yInput[i];

   double halfStep = hstep * 0.5; 

   // Do two half steps
   
   static_cast<T_Stepper*>(this)->T_Stepper::StepWithoutErrorEst (yInitial,  dydx,   halfStep, yMiddle);
   this->RightHandSide(yMiddle, dydxMid);    
   static_cast<T_Stepper*>(this)->T_Stepper::StepWithoutErrorEst (yMiddle, dydxMid, halfStep, yOutput); 

   // Store midpoint, chord calculation

   fMidPoint = ThreeVector( yMiddle[0],  yMiddle[1],  yMiddle[2]); 

   // Do a full Step
   //            static_cast<T_Stepper*>(this)->T_Stepper::StepWithoutErrorEst (yInitial, dydx, hstep, yOneStep);
   static_cast<T_Stepper*>(this)->T_Stepper::StepWithoutErrorEst (yInitial, dydx, hstep, yOneStep);
   for(unsigned int i=0;i<Nvar;i++) {
      yError [i] = yOutput[i] - yOneStep[i] ;
      yOutput[i] += yError[i]* static_cast<T_Stepper*>(this)->T_Stepper::IntegratorCorrection();  
      // T_Stepper::IntegratorCorrection ;
      // Provides accuracy increased by 1 order via the 
      // Richardson Extrapolation  
   }
   
   fInitialPoint = ThreeVector( yInitial[0], yInitial[1], yInitial[2]); 
   fFinalPoint   = ThreeVector( yOutput[0],  yOutput[1],  yOutput[2]); 
   
   return ;
 }


// #ifdef OPT_CHORD_FUNCTIONALITY
template<class T_Stepper, class T_Equation, unsigned int Nvar>
double
TMagErrorStepper<T_Stepper, T_Equation, Nvar>::DistChord() const 
{
            // Estimate the maximum distance from the curve to the chord
            //
            //  We estimate this using the distance of the midpoint to 
            //  chord (the line between 
            // 
            //  Method below is good only for angle deviations < 2 pi, 
            //   This restriction should not a problem for the Runge cutta methods, 
            //   which generally cannot integrate accurately for large angle deviations.
            double distLine, distChord; 

            if (fInitialPoint != fFinalPoint) {
                distLine= GULineSection::Distline( fMidPoint, fInitialPoint, fFinalPoint );
                // This is a class method that gives distance of Mid 
                //  from the Chord between the Initial and Final points.

                distChord = distLine;
            }else{
                distChord = (fMidPoint-fInitialPoint).Mag();
            }

            return distChord;
}
//#endif


#endif  /* TMagErrorStepper_HH */
