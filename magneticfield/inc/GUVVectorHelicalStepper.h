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

#ifndef GUVVectorHelicalStepper_h
#define GUVVectorHelicalStepper_h  1

// #include <CLHEP/Units/PhysicalConstants.h>

// #include "G4Types.hh"

#include "GUVVectorIntegrationStepper.h"
#include "GUVVectorEquationOfMotion.h"
// #include "TMagFieldEquation.h"

// #include "ThreeVector.h"
#include <base/Vector3D.h> 
// typedef vecgeom::Vector3D<double>  ThreeVector; 

class GUVVectorHelicalStepper : public GUVVectorIntegrationStepper
{
  public:  // with description

    typedef Vc::Vector<double> Double_v;
    typedef vecgeom::Vector3D<Double_v>  ThreeVectorSimd; 

    GUVVectorHelicalStepper(GUVVectorEquationOfMotion *EqRhs, // OR TMagFieldEquation *EqRhs,
                            unsigned int order              );
    virtual ~GUVVectorHelicalStepper();
  
    virtual void StepWithErrorEstimate( const Double_v y[],  // virtual for ExactHelix 
                                        const Double_v dydx[],
                                              double h,
                                              Double_v yout[],
                                              Double_v yerr[] );
      // The stepper for the Runge Kutta integration.
      // The stepsize is fixed, equal to h.
      // Integrates ODE starting values y[0 to 6]
      // Outputs yout[] and its estimated error yerr[].
  
    virtual  void StepWithoutErrorEstimate( const Double_v y[],
                                                  ThreeVectorSimd Bfld,
                                                  double  h,
                                                  Double_v yout[]     ) = 0;
      // Performs a 'dump' Step without error calculation.
  
    Double_v DistChord()const ;
      // Estimate maximum distance of curved solution and chord ... 

    virtual void InitializeCharge(double particleCharge)
    {
       fParticleCharge= particleCharge;
       
       // Pass it along as expected 
       GUVVectorIntegrationStepper::InitializeCharge(particleCharge);
    }
       //  GetEquationOfMotion()->InitializeCharge(particleCharge); }
  protected:  // with description

    inline void LinearStep( const Double_v  yIn[],
                                  double  h,
                                  Double_v  yHelix[]) const;
      // A linear Step in regions without magnetic field.

     void AdvanceHelix( const Double_v  yIn[],
                              ThreeVectorSimd Bfld,
                              double  h,
                              Double_v  yHelix[],
                              Double_v yHelix2[]=0);    // output 
      // A first order Step along a helix inside the field.

    inline void MagFieldEvaluate( const Double_v y[], ThreeVectorSimd& Bfield );
      // Evaluate the field at a certain point.

  
    inline Double_v GetInverseCurve( const Double_v Momentum, const Double_v Bmag );
      // Evaluate Inverse of Curvature of Track

      // Store and use the parameters of track : 
      // Radius of curve, Stepping angle, Radius of projected helix
    inline void SetAngCurve(const Double_v Ang);
    inline Double_v GetAngCurve()const;

    inline void SetCurve(const Double_v Curve);
    inline Double_v GetCurve()const;

    inline void SetRadHelix(const Double_v Rad);
    inline Double_v GetRadHelix()const;


  protected:  // without description

    // void MagFieldEvaluate( const double y[], double B[] )   
    //  { GetEquationOfMotion()->  GetFieldValue(y, B); }

  private:

    GUVVectorHelicalStepper(const GUVVectorHelicalStepper&);
    GUVVectorHelicalStepper& operator=(const GUVVectorHelicalStepper&);
      // Private copy constructor and assignment operator.
 
    static const double fUnitConstant;   //  As in TMagFieldEquation.h/cc where it is not used.
  private:
   
      // TMagFieldEquation*
      GUVVectorEquationOfMotion*  fPtrMagEqOfMot;

    // Data stored in order to find the chord.
      Double_v fAngCurve;
      Double_v frCurve;
      Double_v frHelix;
    // Data stored in order to find the chord.
      ThreeVectorSimd yInitial, yMidPoint, yFinal;
       
      double fParticleCharge;
};



inline 
void GUVVectorHelicalStepper::LinearStep( const Vc::Vector<double> yIn[],
                                           double  h,
                                           Vc::Vector<double> yLinear[]) const
{

  typedef Vc::Vector<double> Double_v;
  Double_v  momentum_val = Vc::sqrt(yIn[3]*yIn[3] + yIn[4]*yIn[4] + yIn[5]*yIn[5]) ;
  Double_v  inv_momentum = 1.0 / momentum_val ;
  Double_v  yDir[3];
  // double  h_div_momentum = 1.0 / momentum_val ;

  for( int i = 0; i < 3; i++ ) {
    yDir[i]      = inv_momentum * yIn[i+3];
    yLinear[i]   = yIn[i] + h * yDir[i];
    // yLinear[i]   = yIn[i] + h_div_momentum * yIn[i+3];
    yLinear[i+3] = yIn[i+3];
  }
}

inline 
void GUVVectorHelicalStepper::MagFieldEvaluate(const Vc::Vector<double> y[],
                                                     vecgeom::Vector3D<Vc::Vector<double> > &Bfield )   
{
  Vc::Vector<double> B[3];
  GetEquationOfMotion()->  GetFieldValue(y, B);

  typedef vecgeom::Vector3D<Vc::Vector<double> > ThreeVectorSimd;
  Bfield= ThreeVectorSimd( B[0], B[1], B[2] );
}

inline 
Vc::Vector<double> GUVVectorHelicalStepper::GetInverseCurve(const Vc::Vector<double> Momentum,
                                                            const Vc::Vector<double> Bmag    )   
{
   // define EquationType = TMagFieldEquation<>;
   Vc::Vector<double>  inv_momentum = 1.0 / Momentum ;
   // double particleCharge
   //    = (dynamic_cast<EquationType*>(fPtrMagEqOfMot))->GetParticleCharge(); 
   //     = fPtrMagEqOfMot->FCof() / (CLHEP::eplus*CLHEP::c_light); 
   Vc::Vector<double> fCoefficient = -fUnitConstant * fParticleCharge *inv_momentum;
 
  return  fCoefficient*Bmag;
}


inline 
void GUVVectorHelicalStepper:: SetAngCurve(const Vc::Vector<double> Ang)
{                                                
   fAngCurve=Ang; 
}

inline 
Vc::Vector<double> GUVVectorHelicalStepper:: GetAngCurve() const 
{
  return fAngCurve;  
}

inline 
void GUVVectorHelicalStepper:: SetCurve(const Vc::Vector<double> Curve)
{
  frCurve=Curve;  
}

inline 
Vc::Vector<double> GUVVectorHelicalStepper:: GetCurve() const 
{
  return frCurve;  
}

inline 
void GUVVectorHelicalStepper:: SetRadHelix(const Vc::Vector<double> Rad)
{
  frHelix=Rad;  
}

inline 
Vc::Vector<double> GUVVectorHelicalStepper:: GetRadHelix() const 
{
return frHelix;  

}


#endif  /* GUVVectorHelicalStepper_h */



