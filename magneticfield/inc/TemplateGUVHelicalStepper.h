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

#ifndef TemplateGUVHelicalStepper_h
#define TemplateGUVHelicalStepper_h  1


#include "TemplateGUVIntegrationStepper.h"
#include "TemplateGUVEquationOfMotion.h"

#include <base/Vector3D.h> 

#include "Units.h"     
using fieldUnits::meter;  //  Update to GeantV units ASAP
using fieldUnits::GeV;
using fieldUnits::tesla;

#include "Constants.h"
using Constants::pi;
using Constants::twopi;


template <class Backend>
class TemplateGUVHelicalStepper : public TemplateGUVIntegrationStepper<Backend>
{
  public:  // with description

    typedef typename Backend::precision_v Double_v;
    typedef vecgeom::Vector3D<Double_v>  ThreeVectorSimd; 

    TemplateGUVHelicalStepper(TemplateGUVEquationOfMotion<Backend> *EqRhs, // OR TMagFieldEquation *EqRhs,
                                                        unsigned int order);
    virtual ~TemplateGUVHelicalStepper();
  
    virtual void StepWithErrorEstimate( const Double_v y[],  // virtual for ExactHelix 
                                        const Double_v dydx[],
                                              Double_v h,
                                              Double_v yout[],
                                              Double_v yerr[] );
    // The stepper for the Runge Kutta integration.
    // The stepsize is fixed, equal to h.
    // Integrates ODE starting values y[0 to 6]
    // Outputs yout[] and its estimated error yerr[].
  
    virtual void StepWithoutErrorEstimate( const Double_v y[],
                                                 ThreeVectorSimd Bfld,
                                                 Double_v  h,
                                                 Double_v yout[]     ) = 0;
    // Performs a 'dump' Step without error calculation.
  
    Double_v DistChord()const ;
    // Estimate maximum distance of curved solution and chord ... 

    virtual void InitializeCharge(Double_v particleCharge)
    {
       fParticleCharge= particleCharge;
       
       // Pass it along as expected 
       TemplateGUVIntegrationStepper<Backend>::InitializeCharge(particleCharge);
    }
    //  GetEquationOfMotion()->InitializeCharge(particleCharge); }
  protected:  // with description

    inline void LinearStep( const Double_v  yIn[],
                                  Double_v  h,
                                  Double_v  yHelix[]) const;
    // A linear Step in regions without magnetic field.

     void AdvanceHelix( const Double_v  yIn[],
                              ThreeVectorSimd Bfld,
                              Double_v  h,
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


  private:

    TemplateGUVHelicalStepper(const TemplateGUVHelicalStepper<Backend>&);
    TemplateGUVHelicalStepper& operator=(const TemplateGUVHelicalStepper<Backend>&);
      // Private copy constructor and assignment operator.
 
    static const double fUnitConstant;   //  As in TMagFieldEquation.h/cc where it is not used.
  private:
   
    TemplateGUVEquationOfMotion<Backend>*  fPtrMagEqOfMot;

    // Data stored in order to find the chord.
    Double_v fAngCurve;
    Double_v frCurve;
    Double_v frHelix;
    // Data stored in order to find the chord.
    ThreeVectorSimd yInitial, yMidPoint, yFinal;
       
    Double_v fParticleCharge;
};


template <class Backend>
const double TemplateGUVHelicalStepper<Backend>::fUnitConstant = 0.299792458*(GeV/(tesla*meter));

template <class Backend>
inline 
void TemplateGUVHelicalStepper<Backend>::
  LinearStep( const typename Backend::precision_v yIn[],
                    typename Backend::precision_v h,
                    typename Backend::precision_v yLinear[]) const
{

  typedef typename Backend::precision_v Double_v;
  Double_v momentum_val = vecgeom::Sqrt(yIn[3]*yIn[3] + yIn[4]*yIn[4] + yIn[5]*yIn[5]) ;
  Double_v inv_momentum = 1.0 / momentum_val ;
  Double_v yDir[3];

  for( int i = 0; i < 3; i++ ) {
    yDir[i]      = inv_momentum * yIn[i+3];
    yLinear[i]   = yIn[i] + h * yDir[i];
    yLinear[i+3] = yIn[i+3];
  }
}


template <class Backend>
inline 
void TemplateGUVHelicalStepper<Backend>::
  MagFieldEvaluate(const typename Backend::precision_v y[],
                         vecgeom::Vector3D<typename Backend::precision_v> &Bfield )   
{
  typedef typename Backend::precision_v Double_v;
  Double_v B[3];

  // TemplateGUVEquationOfMotion<Backend> *testEoM = this->GetEquationOfMotion();

  this->GetEquationOfMotion()->  GetFieldValue(y, B);

  typedef vecgeom::Vector3D<Double_v> ThreeVectorSimd;
  Bfield= ThreeVectorSimd( B[0], B[1], B[2] );
}


template <class Backend>
inline 
typename Backend::precision_v 
TemplateGUVHelicalStepper<Backend>::
  GetInverseCurve(const typename Backend::precision_v Momentum,
                  const typename Backend::precision_v Bmag    )   
{
  // define EquationType = TMagFieldEquation<>;
  typename Backend::precision_v inv_momentum = 1.0 / Momentum ;
  // double particleCharge
  //    = (dynamic_cast<EquationType*>(fPtrMagEqOfMot))->GetParticleCharge(); 
  //     = fPtrMagEqOfMot->FCof() / (CLHEP::eplus*CLHEP::c_light); 
  typename Backend::precision_v fCoefficient = -fUnitConstant * fParticleCharge *inv_momentum;
  return  fCoefficient*Bmag;
}


template <class Backend>
inline 
void TemplateGUVHelicalStepper<Backend>::
  SetAngCurve(const typename Backend::precision_v Ang)
{                                                
  fAngCurve=Ang; 
}

template <class Backend>
inline 
typename Backend::precision_v 
TemplateGUVHelicalStepper<Backend>:: 
  GetAngCurve() const 
{
  return fAngCurve;  
}

template <class Backend>
inline 
void TemplateGUVHelicalStepper<Backend>:: 
  SetCurve(const typename Backend::precision_v Curve)
{
  frCurve=Curve;  
}

template <class Backend>
inline 
typename Backend::precision_v 
TemplateGUVHelicalStepper<Backend>:: 
  GetCurve() const 
{
  return frCurve;  
}

template <class Backend>
inline 
void TemplateGUVHelicalStepper<Backend>:: 
  SetRadHelix(const typename Backend::precision_v Rad)
{
  frHelix=Rad;  
}

template <class Backend>
inline 
typename Backend::precision_v TemplateGUVHelicalStepper<Backend>:: 
  GetRadHelix() const 
{
return frHelix;  

}

template <class Backend>
TemplateGUVHelicalStepper<Backend>::
  TemplateGUVHelicalStepper(TemplateGUVEquationOfMotion<Backend> *EqRhs,
                                                 unsigned int order    )            
   : TemplateGUVIntegrationStepper<Backend>(EqRhs, order, 6,  6), //integrate over 6 variables only, // state could be 8 - also t, E
     fPtrMagEqOfMot(EqRhs), 
     fAngCurve(0.), 
     frCurve(0.), 
     frHelix(0.), 
     fParticleCharge(0.0)
{
}

template <class Backend>
TemplateGUVHelicalStepper<Backend>::
  ~TemplateGUVHelicalStepper()
{
}

template <class Backend>
void
TemplateGUVHelicalStepper<Backend>::
  AdvanceHelix( const typename Backend::precision_v  yIn[],
    vecgeom::Vector3D<typename Backend::precision_v > Bfld,    
                      typename Backend::precision_v  h,
                      typename Backend::precision_v  yHelix[],
                      typename Backend::precision_v  yHelix2[])
{
  // const G4int    nvar = 6;
 
  // OLD  const double approc_limit = 0.05;
  // OLD  approc_limit = 0.05 gives max.error=x^5/5!=(0.05)^5/5!=2.6*e-9
  // NEW  approc_limit = 0.005 gives max.error=x^5/5!=2.6*e-14

  // const double approc_limit = 0.005;

  typedef typename Backend::precision_v Double_v;
  typedef vecgeom::Vector3D<typename Backend::precision_v > ThreeVectorSimd;
  // typedef typename Backend::bool_v Bool_v;


  Double_v R_Helix;
  Double_v CosT2, SinT2, CosT, SinT;
  ThreeVectorSimd positionMove, endTangent;

  Double_v Bmag = Bfld.Mag2();
  const Double_v *pIn = yIn+3;
  ThreeVectorSimd initVelocity= ThreeVectorSimd( pIn[0], pIn[1], pIn[2]);
  Double_v        velocityVal = initVelocity.Mag2();
  ThreeVectorSimd initTangent = (1.0/velocityVal) * initVelocity;
  
  Double_v R_1 = GetInverseCurve(velocityVal,Bmag);



  //from if statement
  // LinearStep( yIn, h, yHelix );
  
  //from else statement
  ThreeVectorSimd Bnorm = (1.0/Bmag)*Bfld;
  ThreeVectorSimd B_x_P = Bnorm.Cross(initTangent);
  Double_v        B_d_P = Bnorm.Dot(initTangent); // this is the fraction of P parallel to B
  ThreeVectorSimd vpar  = B_d_P * Bnorm;       // the component parallel      to B
  ThreeVectorSimd vperp = initTangent - vpar;  // the component perpendicular to B
  Double_v        B_v_P = vecgeom::Sqrt( 1 - B_d_P * B_d_P); // Fraction of P perp to B
  Double_v        Theta = R_1 * h; // * B_v_P;


/*  //else inside else
  Double_v Theta2 = Theta  * Theta;
  Double_v Theta3 = Theta2 * Theta;
  Double_v Theta4 = Theta2 * Theta2;
  SinT     = Theta - 1.0/6.0 * Theta3;
  CosT     = 1 - 0.5 * Theta2 + 1.0/24.0 * Theta4;*/

  //if inside else
/*  Bool_v trigCond = vecgeom::Abs(Theta) > approc_limit;
  vecgeom::MaskedAssign(trigCond, vecgeom::sin(Theta), &SinT);
  vecgeom::MaskedAssign(trigCond, vecgeom::cos(Theta), &CosT);*/

  SinT = vecgeom::sin(Theta);
  CosT = vecgeom::cos(Theta);

  Double_v R = 1.0 / R_1;

  positionMove  = R * ( SinT * vperp + (1-CosT) * B_x_P) + h * vpar;
  endTangent    = CosT * vperp + SinT * B_x_P + vpar;

  // Store the resulting position and tangent

  // yHelix[0] = yIn[0] + positionMove.x(); 
  // yHelix[1] = yIn[1] + positionMove.y(); 
  // yHelix[2] = yIn[2] + positionMove.z();
  // yHelix[3] = velocityVal * endTangent.x();
  // yHelix[4] = velocityVal * endTangent.y();
  // yHelix[5] = velocityVal * endTangent.z();

  //try auto-vectorization for above 6 statements:
  for (int i = 0; i < 3; ++i)
  {
    yHelix[i]   = yIn[i] + positionMove[i];
    yHelix[i+3] = velocityVal * endTangent[i];
  }

  //calculations if yHelix2 exists
  SinT2     = 2.0 * SinT * CosT;
  CosT2     = 1.0 - 2.0 * SinT * SinT;
  endTangent    = (CosT2 * vperp + SinT2 * B_x_P + vpar);
  positionMove  = R * ( SinT2 * vperp + (1-CosT2) * B_x_P) + h*2 * vpar;


  for (int i = 0; i < 3; ++i)
  {
    yHelix2[i]   = yIn[i] + positionMove[i];
    yHelix2[i+3] = velocityVal * endTangent[i];
  }


  Double_v ptan=velocityVal*B_v_P;

  R_Helix =vecgeom::Abs( ptan/(fUnitConstant  * fParticleCharge*Bmag));
     

  // for too small magnetic fields there is no curvature
  // (include momentum here) FIXME

/*  Bool_v noCurvatureCond = (vecgeom::Abs(R_1) < 1e-10) || (Bmag<1e-12);
  vecgeom::MaskedAssign(noCurvatureCond, 1., &Theta  );
  vecgeom::MaskedAssign(noCurvatureCond, h , &R      );
  vecgeom::MaskedAssign(noCurvatureCond, 0., &R_Helix);*/


  SetAngCurve(vecgeom::Abs(Theta));
  SetCurve(vecgeom::Abs(R));
  SetRadHelix(R_Helix);


}


//
//  Use the midpoint method to get an error estimate and correction
//  modified from G4ClassicalRK4: W.Wander <wwc@mit.edu> 12/09/97
//

template <class Backend>
void TemplateGUVHelicalStepper<Backend>::
  StepWithErrorEstimate( const typename Backend::precision_v yInput[],
                         const typename Backend::precision_v*,
                               typename Backend::precision_v hstep,
                               typename Backend::precision_v yOut[],
                               typename Backend::precision_v yErr[]  )
{  
   const int nvar = 6;

   typename Backend::precision_v  yTemp[7], yIn[7] ;

   typedef vecgeom::Vector3D<typename Backend::precision_v > ThreeVectorSimd;

   ThreeVectorSimd Bfld_initial, Bfld_midpoint;
   
   //  Saving yInput because yInput and yOut can be aliases for same array

   for(unsigned int i=0;i<nvar;i++) { yIn[i]=yInput[i]; }

   double h = hstep * 0.5; 

   MagFieldEvaluate(yIn, Bfld_initial) ;      

   // Do two half steps

   StepWithoutErrorEstimate(yIn,   Bfld_initial,  h, yTemp);
   MagFieldEvaluate(yTemp, Bfld_midpoint) ;     
   StepWithoutErrorEstimate(yTemp, Bfld_midpoint, h, yOut); 

   // Do a full Step

   h = hstep ;
   StepWithoutErrorEstimate(yIn, Bfld_initial, h, yTemp);

   // Error estimation

   for(unsigned int i=0; i<nvar; ++i)
   {
     yErr[i] = yOut[i] - yTemp[i] ;
   }
   
   return;
}


template <class Backend>
typename Backend::precision_v 
TemplateGUVHelicalStepper<Backend>::
  DistChord() const 
{
  // Check whether h/R >  pi  !!
  // Method DistLine is good only for <  pi

  typename Backend::precision_v Ang=GetAngCurve();
  typename Backend::precision_v returnValue;

  vecgeom::MaskedAssign( Ang<=pi, GetRadHelix()*(1-Vc::cos(0.5*Ang)), &returnValue);
  vecgeom::MaskedAssign( Ang>pi && Ang<twopi, GetRadHelix()*(1+Vc::cos(0.5*(twopi-Ang))), &returnValue );
  vecgeom::MaskedAssign( Ang>= twopi, 2*GetRadHelix(), &returnValue );

  return returnValue;

}


#endif  /* GUVVectorHelicalStepper_h */



