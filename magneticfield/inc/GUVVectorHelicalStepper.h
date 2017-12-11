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

// #include "GUVVectorIntegrationStepper.h"
#include "GUVVectorEquationOfMotion.h"
// #include "ScalarFieldEquation.h"

// #include "ThreeVector.h"
#include <base/Vector3D.h> 
#include <Geant/VectorTypes.h>

class GUVVectorHelicalStepper //  : public GUVVectorIntegrationStepper
{
  using Double_v = Geant::Double_v;
  using Bool_v = Geant::MaskD_v;
  using ThreeVectorSimd = vecgeom::Vector3D<Double_v>;

  public:  // with description

    GUVVectorHelicalStepper(GUVVectorEquationOfMotion *EqRhs, // OR ScalarFieldEquation *EqRhs,
                            unsigned int order              );
    virtual ~GUVVectorHelicalStepper();
  
    virtual void StepWithErrorEstimate( const Double_v y[],  // virtual for ExactHelix 
                                        const Double_v dydx[],
                                        const Double_v&  charge,
                                        const Double_v&  h,                                        
                                              Double_v yout[],
                                              Double_v yerr[] );
      // The stepper for the Runge Kutta integration.
      // The stepsize is fixed, equal to h.
      // Integrates ODE starting values y[0 to 6]
      // Outputs yout[] and its estimated error yerr[].
  
    virtual  void StepWithoutErrorEstimate( const Double_v y[],
                                            ThreeVectorSimd Bfld,
                                            const Double_v& charge,
                                            const Double_v& hStep,
                                                  Double_v yout[]     ) = 0;
      // Performs a 'dump' Step without error calculation.
  
    Double_v DistChord()const ;
      // Estimate maximum distance of curved solution and chord ... 

  protected:  // with description

    inline void LinearStep( const Double_v  yIn[],
                            const Double_v& h,
                                  Double_v  yHelix[]) const;
      // A linear Step in regions without magnetic field.

     void AdvanceHelix( const Double_v  yIn[],
                              ThreeVectorSimd Bfld,
                        const Double_v& charge,
                        const Double_v& h,
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
 
    static const double fUnitConstant;   //  As in ScalarFieldEquation.h/cc where it is not used.
  private:
   
      // ScalarFieldEquation*
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
void GUVVectorHelicalStepper::LinearStep( const Double_v yIn[],
                                          const Double_v& h,
                                          Double_v yLinear[]) const
{

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
void GUVVectorHelicalStepper::MagFieldEvaluate(const Double_v y[],
                                               ThreeVectorSimd &Bfield )
{
  Double_v B[3];
  const char *methodName= "GUVVectorHelicalStepper::MagFieldEvaluate";
  
  auto equation= GetABCEquationOfMotion();
  if( equation ) {
     GetABCEquationOfMotion()->  GetFieldValue(y, B);
  } else {
     std::cerr << "ERROR in " << methodName << ": called when equation held is null\n";
  }
  // GetEquationOfMotion()->  GetFieldValue(y, B);  
  // fPtrMagEqOfMot->GetFieldValue(y, B);

  Bfield.Set(B[0], B[1], B[2]);
}

inline 
Geant::Double_v GUVVectorHelicalStepper::GetInverseCurve(const Double_v Momentum,
                                                  const Double_v Bmag    )   
{
   // define EquationType = ScalarFieldEquation<>;
   Double_v  inv_momentum = 1.0 / Momentum ;
   // double particleCharge
   //    = (dynamic_cast<EquationType*>(fPtrMagEqOfMot))->GetParticleCharge(); 
   //     = fPtrMagEqOfMot->FCof() / (CLHEP::eplus*CLHEP::c_light); 
   Double_v coefficient = -fUnitConstant * fParticleCharge *inv_momentum;
 
  return  coefficient*Bmag;
}


inline 
void GUVVectorHelicalStepper:: SetAngCurve(const Double_v Ang)
{                                                
   fAngCurve=Ang; 
}

inline 
Geant::Double_v GUVVectorHelicalStepper:: GetAngCurve() const 
{
  return fAngCurve;  
}

inline 
void GUVVectorHelicalStepper:: SetCurve(const Double_v Curve)
{
  frCurve=Curve;  
}

inline 
Geant::Double_v GUVVectorHelicalStepper:: GetCurve() const 
{
  return frCurve;  
}

inline 
void GUVVectorHelicalStepper:: SetRadHelix(const Double_v Rad)
{
  frHelix=Rad;  
}

inline 
Geant::Double_v GUVVectorHelicalStepper:: GetRadHelix() const 
{
return frHelix;  

}


#endif  /* GUVVectorHelicalStepper_h */
