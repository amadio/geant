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
                               int numStateVariables); // = -1 same? or  unsigned ?    // in G4 =12
           // See explanations of each below - e.g. order => RK order

        GUVIntegrationStepper( const GUVIntegrationStepper& );
           // For use in Clone() method
        
        virtual ~GUVIntegrationStepper();

        // Core methods
        // ---------------------
        virtual void StepWithErrorEstimate( const double y[],
                                            const double dydx[],
                                                  double charge,
                                                  double h,
                                                  double yout[],
                                                  double yerr[]  ) = 0;
        // Integrate typically using Runge Kutta 
        // Input:
        //          y[] = initial derivative
        //       dydx[] = initial derivative        
        //          h   = requested step
        // Output:
        //       yout[] = output values of integration
        //       yerr[] = estimate of integration error

        virtual  double  DistChord(double charge) const = 0; 
        // Estimate the maximum sagital distance (distance of a chord from the true path)
        //  over the last segment integrated.

        // Auxiliary methods
        // ---------------------
        virtual  GUVIntegrationStepper* Clone() const = 0;
        // Create an independent copy of the current object -- including independent 'owned' objects
        
        inline void RightHandSideVIS( const double y[], double charge, double dydx[] );   
        // Utility method to supply the standard Evaluation of the
        // Right Hand side of the associated equation.

        // virtual void ComputeRightHandSide( const double y[], double charge, double dydx[] ); 
        // Must compute the RightHandSide as in the method above
        // Optionally can cache the input y[] and the dydx[] values computed.

        inline unsigned int  GetNumberOfVariables() const;
        
        // Get the number of variables that the stepper will integrate over.

        inline unsigned int  GetNumberOfStateVariables() const;
        // Get the number of variables of state variables (>= above, integration)

        unsigned int IntegratorOrder() const { return fIntegrationOrder; };
        // Returns the order of the integrator
        // i.e. its error behaviour is of the order O(h^order).

        // inline void NormalisePolarizationVector( double vec[12] ); // TODO - add polarisation
        // Simple utility function to (re)normalise 'unit spin' vector.

        inline GUVEquationOfMotion *GetEquationOfMotion() { return  fAbstrEquation; }
        inline const GUVEquationOfMotion *GetEquationOfMotion() const { return  fAbstrEquation; }        
        // As some steppers require access to other methods of Eq_of_Mot
        void SetEquationOfMotion(GUVEquationOfMotion* newEquation); 

//        virtual void InitializeCharge(double particleCharge) {
//               GetEquationOfMotion()->InitializeCharge(particleCharge); }
           // Some steppers may need the value(s) / or status - they can intercept        

    private:

        GUVIntegrationStepper& operator=(const GUVIntegrationStepper&);
        // Private copy constructor and assignment operator.

    private:

        GUVEquationOfMotion *fAbstrEquation;  // For use in calling RightHandSideVIS only
          // Object is typically owned by stepper, but if a separate pointer (TEquation)
          //  exists which points to the same object, it must not be deleted using
          //  this pointer!
        
        const unsigned int fIntegrationOrder; // RK or similar order - if any. Else 0
        const unsigned int fNoIntegrationVariables; // # of Variables in integration
        const unsigned int fNoStateVariables;       // # required for FieldTrack
};

// #include  "GUVIntegrationStepper.icc"
inline
void GUVIntegrationStepper::
RightHandSideVIS( const  double y[], double charge, double dydx[] )
{
   assert ( fAbstrEquation != nullptr ); 
   fAbstrEquation-> RightHandSide(y, charge, dydx);
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
