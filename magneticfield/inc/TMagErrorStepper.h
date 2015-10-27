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
<class T_Stepper, class T_Equation, int Nvar>
class TMagErrorStepper : public GUVIntegrationStepper
{
    public:  // with description
        static constexpr unsigned int NumVarStore = (Nvar > GUIntegrationNms::NumVarBase) ?
                                                     Nvar : GUIntegrationNms::NumVarBase ;
           // std::max( GUIntegrationNms::NumVarBase,  Nvar);

        TMagErrorStepper( T_Equation *EqRhs,
                          unsigned int integrationOrder,   // Make it a template Parameter ??
                          unsigned int numberOfVariables,   // Redundant -- to be DELETED
                          int numStateVariables) // = -1)  // No default -- must ensure order is set
            : GUVIntegrationStepper(
                  EqRhs,
                  integrationOrder,      
                  numberOfVariables,
                  numStateVariables ), // ((numStateVariables>0) ? numStateVariables : NumVarStore) ),
               // ),
               fEquation_Rhs(EqRhs)
        {
            // int nvar = std::max(this->GetNumberOfVariables(), 8);
        }

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

        TMagErrorStepper(const TMagErrorStepper&);
        TMagErrorStepper& operator=(const TMagErrorStepper&);
        // Private copy constructor and assignment operator.

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

// inline
template<class T_Stepper, class T_Equation, int Nvar>
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
            const int maxvar= GetNumberOfStateVariables();

            int i;
            // correction for Richardson Extrapolation.
            //double  correction = 1. / ( (1 << 
              //          static_cast<T_Stepper*>(this)->T_Stepper::IntegratorOrder()) -1 );
            //  Saving yInput because yInput and yOutput can be aliases for same array

            for(i=0;i<N;i++) yInitial[i]=yInput[i];
            yInitial[7]= yInput[7];    // Copy the time in case ... even if not really needed
            yMiddle[7] = yInput[7];  // Copy the time from initial value 
            yOneStep[7] = yInput[7]; // As it contributes to final value of yOutput ?
            // yOutput[7] = yInput[7];  // -> dumb stepper does it too for RK4
            for(i=N;i<maxvar;i++) yOutput[i]=yInput[i];
            // yError[7] = 0.0;         

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
            for(i=0;i<Nvar;i++) {
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

#ifdef OPT_CHORD_FUNCTIONALITY        
        double DistChord() const 
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
#endif

    private:

        TMagErrorStepper(const TMagErrorStepper&);
        TMagErrorStepper& operator=(const TMagErrorStepper&);
        // Private copy constructor and assignment operator.

    private:

        // STATE
        ThreeVector fInitialPoint, fMidPoint, fFinalPoint;
        // Data stored in order to find the chord

        // Dependent Objects, owned --- part of the STATE 
        double yInitial[N<8?8:N];
        double yMiddle[N<8?8:N];
        double dydxMid[N<8?8:N];
        double yOneStep[N<8?8:N];
        // The following arrays are used only for temporary storage
        // they are allocated at the class level only for efficiency -
        // so that calls to new and delete are not made in Stepper().

        T_Equation *fEquation_Rhs;
};


#endif  /* TMagErrorStepper_HH */
