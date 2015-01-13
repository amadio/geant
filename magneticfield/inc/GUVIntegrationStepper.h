//
// class GUVIntegrationStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field
//
//  A Stepper must integrate over                NumberOfVariables elements,
//   and also copy (from input to output) any of NoStateVariables  
//   not included in the NumberOfVariables.  
// 
//  So it is expected that NoStateVariables >= NumberOfVariables

// History:
// - 15.01.97  J. Apostolakis (J.Apostolakis@cern.ch)
// --------------------------------------------------------------------

#ifndef GUVIntegrationStepper_h
#define GUVIntegrationStepper_h

// #include "GUVTypes.h"
#include "GUVEquationOfMotion.h"
// class GUVEquationOfMotion;
 
class GUVIntegrationStepper
{
  public:
        GUVIntegrationStepper();   // DELET
        GUVIntegrationStepper( GUVEquationOfMotion* equation, 
                int              numIntegrationVariables,
                int              numStateVariables=12);
        virtual ~GUVIntegrationStepper();

        virtual  void  Step(  const double y[],
                const double dydx[],
                double h,
                double yout[],
                double yerr[]  ) = 0;
        // The stepper for the Runge Kutta integration.
        // The stepsize is fixed, with the Step size given by h.
        // Integrates ODE starting values y[0 to 6].
        // Outputs yout[] and its estimated error yerr[].

        virtual  double  DistChord() const = 0; 
        // Estimate the maximum distance of a chord from the true path
        // over the segment last integrated.

        virtual void ComputeRightHandSide( const double y[], double charge, double dydx[] ); 
        // Must compute the RightHandSide as in the method below
        // Optionally can cache the input y[] and the dydx[] values computed.

        // inline void NormaliseTangentVector( double vec[6] );  // WRONG - it is Momentum now!!
        // Simple utility function to (re)normalise 'unit velocity' vector.

        // inline void NormalisePolarizationVector( double vec[12] ); // TODO - add polarisation
        // Simple utility function to (re)normalise 'unit spin' vector.

        inline void RightHandSide( const double y[], double charge, double dydx[] );   
        // Utility method to supply the standard Evaluation of the
        // Right Hand side of the associated equation.


        inline int  GetNumberOfVariables() const;
        // Get the number of variables that the stepper will integrate over.

        // void   SetNumberOfVariables(int newNo);  // Dangerous & obsolete ...

        inline int  GetNumberOfStateVariables() const;
        // Get the number of variables of state variables (>= above, integration)

        virtual int IntegratorOrder() const = 0;
        // Returns the order of the integrator
        // i.e. its error behaviour is of the order O(h^order).

        inline GUVEquationOfMotion *GetEquationOfMotion(); 
        // As some steppers (eg RKG3) require other methods of Eq_Rhs
        // this function allows for access to them.
        inline void SetEquationOfMotion(GUVEquationOfMotion* newEquation); 

    private:

        GUVIntegrationStepper(const GUVIntegrationStepper&);
        GUVIntegrationStepper& operator=(const GUVIntegrationStepper&);
        // Private copy constructor and assignment operator.

    private:

        GUVEquationOfMotion *fEquation_Rhs;
        const int  fNoIntegrationVariables;  // Number of Variables in integration
        const int  fNoStateVariables;        // Number required for FieldTrack
        // const int  fNumberOfVariables;
};

#include  "GUVIntegrationStepper.icc"

#endif  /* GUVIntegrationStepper */
