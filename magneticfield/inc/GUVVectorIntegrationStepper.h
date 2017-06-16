//
// class GUVVectorIntegrationStepper
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

#ifndef GUVVectorIntegrationStepper_h
#define GUVVectorIntegrationStepper_h

// #include "GUVTypes.h"
#include "GUVVectorEquationOfMotion.h"
// class GUVVectorEquationOfMotion;
 
class GUVVectorIntegrationStepper
{
    typedef typename vecgeom::kVc::precision_v Double_v;

  public:
        // GUVVectorIntegrationStepper();   // DELET
        GUVVectorIntegrationStepper( GUVVectorEquationOfMotion* equation, 
                                     unsigned int IntegrationOrder,
                                     unsigned int numIntegrationVariables,
                                     int          numStateVariables      ); // = -1 same? or  unsigned ?    // in G4 =12
           // See explanations of each below - e.g. order => RK order

        GUVVectorIntegrationStepper( const GUVVectorIntegrationStepper& );
           // For use in Clone() method
        
        virtual ~GUVVectorIntegrationStepper();

        // Core methods
        // ---------------------
        virtual void StepWithErrorEstimate( const Double_v y[],
                                            const Double_v dydx[],
//                                            const Double_v h,
                                                  double   h, 
                                                  Double_v yout[],
                                                  Double_v yerr[]  ) = 0;
        // Integrate typically using Runge Kutta 
        // Input:
        //          y[] = initial derivative
        //       dydx[] = initial derivative        
        //          h   = requested step
        // Output:
        //       yout[] = output values of integration
        //       yerr[] = estimate of integration error

        virtual  Double_v  DistChord() const = 0; 
        // Estimate the maximum sagital distance (distance of a chord from the true path)
        //  over the last segment integrated.

        // Auxiliary methods
        // ---------------------
        virtual  GUVVectorIntegrationStepper* Clone() const = 0;
        // Create an independent copy of the current object -- including independent 'owned' objects
        
        inline void RightHandSideVIS( const Double_v y[], Double_v charge, Double_v dydx[] );   
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

        inline GUVVectorEquationOfMotion *GetEquationOfMotion() { return  fAbstrEquation; }
        inline const GUVVectorEquationOfMotion *GetEquationOfMotion() const { return  fAbstrEquation; }        
        // As some steppers require access to other methods of Eq_of_Mot
        void SetEquationOfMotion(GUVVectorEquationOfMotion* newEquation); 

        virtual void InitializeCharge(double particleCharge) { GetEquationOfMotion()->InitializeCharge(particleCharge); }
           // Some steppers may need the value(s) / or status - they can intercept        

        void InformDone() { GetEquationOfMotion()->InformDone();}
          // InvalidateParameters()

    private:

        GUVVectorIntegrationStepper& operator=(const GUVVectorIntegrationStepper&);
        // Private copy constructor and assignment operator.

    private:

        GUVVectorEquationOfMotion *fAbstrEquation;  // For use in calling RightHandSideVIS only
          // Object is typically owned by stepper, but if a separate pointer (TEquation)
          //  exists which points to the same object, it must not be deleted using
          //  this pointer!
        
        const unsigned int fIntegrationOrder; // RK or similar order - if any. Else 0
        const unsigned int fNoIntegrationVariables; // # of Variables in integration
        const unsigned int fNoStateVariables;       // # required for FieldTrack
};

// #include  "GUVVectorIntegrationStepper.icc"
inline
void GUVVectorIntegrationStepper::
RightHandSideVIS( const  typename vecgeom::kVc::precision_v y[], 
                         typename vecgeom::kVc::precision_v charge,
                         typename vecgeom::kVc::precision_v dydx[] )
{
   fAbstrEquation-> RightHandSide(y, charge, dydx);
}

inline
unsigned int GUVVectorIntegrationStepper::GetNumberOfVariables() const
{
  return fNoIntegrationVariables;
}


inline
unsigned int GUVVectorIntegrationStepper::GetNumberOfStateVariables() const
{
  return fNoStateVariables;
}

#endif  /* GUVVectorIntegrationStepper */
