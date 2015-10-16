//
//  GUExactHelixStepper
//  -------------------
//
//  Based on G4ExactHelixStepper
//
// Adapted from G4ExactHelixStepper 
// - 16.Oct.15  J.Apostolakis   Adapted
// --------------------------------------------------------------------

#ifndef GUExactHelixStepper_h
#define GUExactHelixStepper_h 1

// #include "G4Types.h"
 #include "ThreeVector.h"

// #include "GUVIntegrationStepper.h"
#include "GUVHelicalStepper.h"
#include "TMagFieldEquation.h"

class GUExactHelixStepper : public GUVHelicalStepper
{
  public:  // with description

    GUExactHelixStepper(GUVEquationOfMotion *EqRhs); // TMagFieldEquation *EqRhs);
    ~GUExactHelixStepper();

    void StepWithErrorEstimate( const double y[],
                  const double dydx[],
                        double h,
                        double yout[],
                        double yerr[]  );
      // Step 'integration' for step size 'h'
      // Provides helix starting at y[0 to 6]
      // Outputs yout[] and ZERO estimated error yerr[]=0.
  
    void StepWithoutErrorEstimate( const double y[],
                            ThreeVector   Bfld,
                            double  h,
                            double yout[] );
      // Performs a 'dump' Step without error calculation.
  
    double DistChord() const;
      // Estimate maximum distance of curved solution and chord ... 

  private:

    GUExactHelixStepper(const GUExactHelixStepper&);
    GUExactHelixStepper& operator=(const GUExactHelixStepper&);
      // Private copy constructor and assignment operator.
   
  private:

    ThreeVector    fBfieldValue;
      //  Initial value of field at last step
    // TMagFieldEquation*
    GUVEquationOfMotion* fPtrMagEqOfMot;
};



#endif  /* GUExactHelixStepper_h */
