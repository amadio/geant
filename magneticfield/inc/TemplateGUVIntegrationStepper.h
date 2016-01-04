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

#ifndef TemplateGUVIntegrationStepper_h
#define TemplateGUVIntegrationStepper_h

#include "TemplateGUVEquationOfMotion.h"

// #define DEBUGAnanya

template <class Backend> 
class TemplateGUVIntegrationStepper
{
    typedef typename Backend::precision_v Double_v;

  public:
        // TemplateGUVIntegrationStepper();   // DELET
        TemplateGUVIntegrationStepper(TemplateGUVEquationOfMotion<Backend>* equation, 
                                                unsigned int IntegrationOrder,
                                                unsigned int numIntegrationVariables,
                                                int          numStateVariables      ); // = -1 same? or  unsigned ?    // in G4 =12
           // See explanations of each below - e.g. order => RK order

        TemplateGUVIntegrationStepper( const TemplateGUVIntegrationStepper& );
           // For use in Clone() method
        
        virtual ~TemplateGUVIntegrationStepper();

        // Core methods
        // ---------------------
        virtual void StepWithErrorEstimate( const Double_v y[],
                                            const Double_v dydx[],
                                            const Double_v charge,
                                                  Double_v h,
                                                  Double_v yout[],
                                                  Double_v yerr[]  ) = 0;
        // Integrate typically using Runge Kutta 
        // Input:
        //          y[] = initial derivative
        //       dydx[] = initial derivative
        //       charge = charge
        //          h   = requested step
        // Output:
        //       yout[] = output values of integration
        //       yerr[] = estimate of integration error

        virtual  Double_v  DistChord() const = 0; 
        // Estimate the maximum sagital distance (distance of a chord from the true path)
        //  over the last segment integrated.

        // Auxiliary methods
        // ---------------------
        virtual  TemplateGUVIntegrationStepper* Clone() const = 0;
        // Create an independent copy of the current object -- including independent 'owned' objects
        
        inline void RightHandSideVIS( const Double_v y[], Double_v charge, Double_v dydx[] );   

        // inline void RightHandSideVIS( const Double_v y[], Double_v charge, Double_v dydx[] ); 

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

        inline TemplateGUVEquationOfMotion<Backend> *GetEquationOfMotion() { return  fAbstrEquation; }
        inline const TemplateGUVEquationOfMotion<Backend> *GetEquationOfMotion() const { return  fAbstrEquation; }        
        // As some steppers require access to other methods of Eq_of_Mot
        void SetEquationOfMotion(TemplateGUVEquationOfMotion<Backend>* newEquation); 

// virtual void InitializeCharge(double particleCharge) { GetEquationOfMotion()->InitializeCharge(particleCharge); }
           // Some steppers may need the value(s) / or status - they can intercept        

    private:

        TemplateGUVIntegrationStepper& operator=(const TemplateGUVIntegrationStepper&);
        // Private copy constructor and assignment operator.

    private:

        TemplateGUVEquationOfMotion<Backend> *fAbstrEquation;  // For use in calling RightHandSideVIS only
          // Object is typically owned by stepper, but if a separate pointer (TEquation)
          //  exists which points to the same object, it must not be deleted using
          //  this pointer!
        
        const unsigned int fIntegrationOrder; // RK or similar order - if any. Else 0
        const unsigned int fNoIntegrationVariables; // # of Variables in integration
        const unsigned int fNoStateVariables;       // # required for FieldTrack
};

// #include  "TemplateGUVIntegrationStepper.icc"
template <class Backend> 
inline
void TemplateGUVIntegrationStepper<Backend>::
  RightHandSideVIS( const  typename Backend::precision_v y[], 
                           typename Backend::precision_v charge,
                           typename Backend::precision_v dydx[] )
{
   fAbstrEquation-> RightHandSide(y, charge, dydx);
/*   
#ifdef DEBUGAnanya
   std::cout<<"\n----y to RightHandSideVIS is: "<<y[3]<<std::endl;
   #endif 
*/
}

template <class Backend> 
inline
unsigned int TemplateGUVIntegrationStepper<Backend>::GetNumberOfVariables() const
{
  return fNoIntegrationVariables;
}

template <class Backend> 
inline
unsigned int TemplateGUVIntegrationStepper<Backend>::GetNumberOfStateVariables() const
{
  return fNoStateVariables;
}


template <class Backend> 
TemplateGUVIntegrationStepper<Backend>::
  TemplateGUVIntegrationStepper(TemplateGUVEquationOfMotion<Backend>* equation,
                                            unsigned int num_integration_vars,
                                            unsigned int integrationOrder,
                                                     int num_state_vars       )
  : fAbstrEquation(equation),
    fIntegrationOrder(integrationOrder),
    fNoIntegrationVariables(num_integration_vars),
    fNoStateVariables(num_state_vars > 0 ? num_state_vars : num_integration_vars)
{
  #if 0
    std::cout<<"\n----Entered constructor of TemplateGUVIntegrationStepper"<<std::endl;
    std::cout<<"----Equation is: "<<equation->idxTime<<std::endl;
    std::cout<<"----num_state_vars is: "<<integrationOrder<<std::endl;
  #endif
}

template <class Backend> 
TemplateGUVIntegrationStepper<Backend>::~TemplateGUVIntegrationStepper()
{
  #if 0
    std::cout<<"----IntegrationStepper destructor"<<std::endl;
  #endif 
}

// This allows the method to cache the value etc - Not needed for now
// void TemplateGUVIntegrationStepper<Backend>::ComputeRightHandSide( const double y[], /*double charge,*/ double dydx[] ) 
// {
//    this->RightHandSide( y, /*charge,*/ dydx );
// }

template <class Backend> 
void TemplateGUVIntegrationStepper<Backend>::
  SetEquationOfMotion(TemplateGUVEquationOfMotion<Backend>* newEquation)
{
  if( newEquation != 0 )
  {
    fAbstrEquation= newEquation;
  }
}


#endif  /* TemplateGUVIntegrationStepper */
