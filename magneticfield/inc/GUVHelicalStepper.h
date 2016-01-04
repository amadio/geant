//
// Based on G4MagHelicalStepper
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field
//
// It is used for a set of steppers which use the helix as a sort of
// 'first order' solution.
//   - Most obtain an error by breaking up the step in two
//   - G4ExactHelicalStepper does not provide an error estimate
//
// History:
// - 16.10.15  J.Apostolakis   Created for testing other steppers
// --------------------------------------------------------------------

#ifndef GUVHelicalStepper_h
#define GUVHelicalStepper_h  1

// #include <CLHEP/Units/PhysicalConstants.h>

// #include "G4Types.hh"

#include "GUVIntegrationStepper.h"
#include "GUVEquationOfMotion.h"
// #include "TMagFieldEquation.h"

// #include "ThreeVector.h"
#include <base/Vector3D.h> 
typedef vecgeom::Vector3D<double>  ThreeVector; 

class GUVHelicalStepper : public GUVIntegrationStepper
{
  public:  // with description

    GUVHelicalStepper(GUVEquationOfMotion *EqRhs, // OR TMagFieldEquation *EqRhs,
                      unsigned int order);
    virtual ~GUVHelicalStepper();
  
    virtual void StepWithErrorEstimate( const double y[],  // virtual for ExactHelix 
                  const double dydx[],
                        double h,
                        double yout[],
                        double yerr[]  );
      // The stepper for the Runge Kutta integration.
      // The stepsize is fixed, equal to h.
      // Integrates ODE starting values y[0 to 6]
      // Outputs yout[] and its estimated error yerr[].
  
    virtual  void StepWithoutErrorEstimate( const double y[],
                               ThreeVector   Bfld,
                               double  h,
                               double yout[] ) = 0;
      // Performs a 'dump' Step without error calculation.
  
    double DistChord()const ;
      // Estimate maximum distance of curved solution and chord ... 

    virtual void InitializeCharge(double particleCharge)
    {
       fParticleCharge= particleCharge;
       
       // Pass it along as expected 
       GUVIntegrationStepper::InitializeCharge(particleCharge);
    }
       //  GetEquationOfMotion()->InitializeCharge(particleCharge); }
  protected:  // with description

    inline void LinearStep( const double  yIn[],
                                  double  h,
                                  double  yHelix[]) const;
      // A linear Step in regions without magnetic field.

     void AdvanceHelix( const double  yIn[],
                             ThreeVector   Bfld,
                             double  h,
                        double  yHelix[],double yHelix2[]=0);    // output 
      // A first order Step along a helix inside the field.

    inline void MagFieldEvaluate( const double y[], ThreeVector& Bfield );
      // Evaluate the field at a certain point.

  
   inline double GetInverseCurve( const double Momentum, const double Bmag );
      // Evaluate Inverse of Curvature of Track

      // Store and use the parameters of track : 
      // Radius of curve, Stepping angle, Radius of projected helix
   inline void SetAngCurve(const double Ang);
   inline double GetAngCurve()const;

   inline void SetCurve(const double Curve);
   inline double GetCurve()const;

   inline void SetRadHelix(const double Rad);
   inline double GetRadHelix()const;


  protected:  // without description

    // void MagFieldEvaluate( const double y[], double B[] )   
    //  { GetEquationOfMotion()->  GetFieldValue(y, B); }

  private:

    GUVHelicalStepper(const GUVHelicalStepper&);
    GUVHelicalStepper& operator=(const GUVHelicalStepper&);
      // Private copy constructor and assignment operator.
 
    static const double fUnitConstant;   //  As in TMagFieldEquation.h/cc where it is not used.
  private:
   
      // TMagFieldEquation*
      GUVEquationOfMotion*  fPtrMagEqOfMot;

    // Data stored in order to find the chord.
      double fAngCurve;
      double frCurve;
      double frHelix;
    // Data stored in order to find the chord.
      ThreeVector yInitial, yMidPoint, yFinal;
       
      double fParticleCharge;
};

// #include  "GUVHelicalStepper.icc"

inline void
GUVHelicalStepper::LinearStep( const double  yIn[],
                                       double  h,
                                       double  yLinear[]) const
{
  double  momentum_val = std::sqrt(yIn[3]*yIn[3] + yIn[4]*yIn[4] + yIn[5]*yIn[5]) ;
  double  inv_momentum = 1.0 / momentum_val ;
  double  yDir[3];
  // double  h_div_momentum = 1.0 / momentum_val ;

  for( int i = 0; i < 3; i++ ) {
    yDir[i]   = inv_momentum * yIn[i+3];
    yLinear[i]   = yIn[i] + h * yDir[i];
    // yLinear[i]   = yIn[i] + h_div_momentum * yIn[i+3];
    yLinear[i+3] = yIn[i+3];
  }
}

inline void
GUVHelicalStepper::MagFieldEvaluate(const double y[],
                                      ThreeVector& Bfield )   
{
  double B[3];
  GetEquationOfMotion()->  GetFieldValue(y, B);
  Bfield= ThreeVector( B[0], B[1], B[2] );
}

inline double
GUVHelicalStepper::GetInverseCurve( const double Momentum,
                                    const double Bmag)   
{
   // define EquationType = TMagFieldEquation<>;
   double  inv_momentum = 1.0 / Momentum ;
   // double particleCharge
   //    = (dynamic_cast<EquationType*>(fPtrMagEqOfMot))->GetParticleCharge(); 
   //     = fPtrMagEqOfMot->FCof() / (CLHEP::eplus*CLHEP::c_light); 
   double fCoefficient = -fUnitConstant * fParticleCharge *inv_momentum;
 
  return  fCoefficient*Bmag;
}
inline void GUVHelicalStepper:: SetAngCurve(const double Ang)
{                                                
   fAngCurve=Ang; 
}

inline double GUVHelicalStepper:: GetAngCurve() const 
{
  return fAngCurve;  
}

inline void GUVHelicalStepper:: SetCurve(const double Curve)
{
  frCurve=Curve;  
}

inline double GUVHelicalStepper:: GetCurve() const 
{
  return frCurve;  
}

inline void GUVHelicalStepper:: SetRadHelix(const double Rad)
{
  frHelix=Rad;  
}

inline double GUVHelicalStepper:: GetRadHelix() const 
{
return frHelix;  

}

#endif  /* GUVHelicalStepper_h */



