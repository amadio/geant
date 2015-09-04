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
// [ So the following must hold: NoStateVariables >= NumberOfVariables ] 
//
//  The integration order is property of convergence of deviation / error,
//   and is meant to be used for (or correspond to) the order of RK method.
//
// First version/origin:
// - Jan-Mar 2015 Created by J. Apostolakis (J.Apostolakis@cern.ch)
//                Derived from my G4MagIntegrationStepper class 
// --------------------------------------------------------------------

#ifndef GUVIntegrationStepper_h
#define GUVIntegrationStepper_h

// #include "GUVTypes.h"
#include "GUVEquationOfMotion.h"
// class GUVEquationOfMotion;
 
class GUVIntegrationStepper
{
  public:
        // GUVIntegrationStepper();   // DELET
        GUVIntegrationStepper( GUVEquationOfMotion* equation, 
                               unsigned int IntegrationOrder,
                               unsigned int numIntegrationVariables,
                               unsigned int numStateVariables); // =12);
           // See explanations of each below - e.g. order => RK order

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


        inline unsigned int  GetNumberOfVariables() const;
        
        // Get the number of variables that the stepper will integrate over.

        // void   SetNumberOfVariables(int newNo);  // Dangerous & obsolete ...

        inline unsigned int  GetNumberOfStateVariables() const;
        // Get the number of variables of state variables (>= above, integration)

        unsigned int IntegratorOrder() { return fIntegrationOrder; };
        // Returns the order of the integrator
        // i.e. its error behaviour is of the order O(h^order).

        inline GUVEquationOfMotion *GetEquationOfMotion() { return  fEquation_Rhs; }
        // As some steppers require access to other methods of Eq_of_Mot
        void SetEquationOfMotion(GUVEquationOfMotion* newEquation); 

    private:

        GUVIntegrationStepper(const GUVIntegrationStepper&);
        GUVIntegrationStepper& operator=(const GUVIntegrationStepper&);
        // Private copy constructor and assignment operator.

    private:

        GUVEquationOfMotion *fEquation_Rhs;
        const unsigned int fIntegrationOrder; // RK or similar order - if any. Else 0
        const unsigned int fNoIntegrationVariables; // # of Variables in integration
        const unsigned int fNoStateVariables;       // # required for FieldTrack
};

// #include  "GUVIntegrationStepper.icc"
inline
void GUVIntegrationStepper::
RightHandSide( const  double y[], double charge, double dydx[] )
{
   fEquation_Rhs-> RightHandSide(y, charge, dydx);
}

inline void
GUVIntegrationStepper::SetEquationOfMotion(GUVEquationOfMotion* newEquation)
{
  if( newEquation != 0 )
  {
    fEquation_Rhs= newEquation;
  }
}

inline
unsigned int GUVIntegrationStepper::GetNumberOfVariables() const
{
  return fNoIntegrationVariables;
}


inline
unsigned int GUVIntegrationStepper::GetNumberOfStateVariables() const
{
  return fNoStateVariables;
}

#endif  /* GUVIntegrationStepper */
