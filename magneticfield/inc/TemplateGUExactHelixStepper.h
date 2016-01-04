//
//  GUExactHelixStepper
//  -------------------
//
//  Based on G4ExactHelixStepper
//
// Adapted from G4ExactHelixStepper 
// - 16.Oct.15  J.Apostolakis   Adapted
// --------------------------------------------------------------------

#ifndef TemplateGUExactHelixStepper_h
#define TemplateGUExactHelixStepper_h 1

//  #include "ThreeVector.h"
#include "base/Vector3D.h"

// #include "GUVIntegrationStepper.h"
#include "TemplateGUVHelicalStepper.h"
#include "TemplateTMagFieldEquation.h"

template <class Backend>
class TemplateGUExactHelixStepper : public TemplateGUVHelicalStepper<Backend>
{
  public:  // with description
    typedef typename Backend::precision_v  Double_v;

    TemplateGUExactHelixStepper(TemplateGUVEquationOfMotion<Backend> *EqRhs); // TMagFieldEquation *EqRhs);
    ~TemplateGUExactHelixStepper();

    void StepWithErrorEstimate( const Double_v y[],
                                const Double_v dydx[],
                                      Double_v h,
                                      Double_v yout[],
                                      Double_v yerr[] );
      // Step 'integration' for step size 'h'
      // Provides helix starting at y[0 to 6]
      // Outputs yout[] and ZERO estimated error yerr[]=0.
  
    void StepWithoutErrorEstimate( const Double_v y[],
                                         vecgeom::Vector3D<Double_v> Bfld,
                                         Double_v  h,
                                         Double_v yout[]                 );
      // Performs a 'dump' Step without error calculation.
  
    Double_v DistChord() const;
      // Estimate maximum distance of curved solution and chord ... 

  private:

    TemplateGUExactHelixStepper(const TemplateGUExactHelixStepper&);
    TemplateGUExactHelixStepper& operator=(const TemplateGUExactHelixStepper&);
      // Private copy constructor and assignment operator.
   
  private:

    // ThreeVector    fBfieldValue;
    vecgeom::Vector3D<Double_v>   fBfieldValue;
      //  Initial value of field at last step

  // TMagFieldEquation*
    TemplateGUVEquationOfMotion<Backend>* fPtrMagEqOfMot;
};


#include "Constants.h"
using Constants::pi;
using Constants::twopi;

template <class Backend>
TemplateGUExactHelixStepper<Backend>
  ::TemplateGUExactHelixStepper(TemplateGUVEquationOfMotion<Backend>* EqRhs) // TMagFieldEquation
   : TemplateGUVHelicalStepper<Backend>(EqRhs, 1),  // "Order" = 1 - not really applicable
     fBfieldValue(DBL_MAX), // , DBL_MAX, DBL_MAX),
     fPtrMagEqOfMot(EqRhs)
{
  ;
}

template <class Backend>
TemplateGUExactHelixStepper<Backend>
  ::~TemplateGUExactHelixStepper() {} 

template <class Backend>
void
TemplateGUExactHelixStepper<Backend>
  ::StepWithErrorEstimate( const typename Backend::precision_v yInput[],
                           const typename Backend::precision_v*,
                                 typename Backend::precision_v hstep,
                                 typename Backend::precision_v yOut[],
                                 typename Backend::precision_v yErr[]  )
{  
   const unsigned int nvar = 6;

   vecgeom::Vector3D<typename Backend::precision_v> Bfld_value;

   MagFieldEvaluate(yInput, Bfld_value);
   // std::cout << " Exact Helix: B-field:  Bx = " << Bfld_value[0]
   //           << " By= " << Bfld_value[1] << " Bz= " << Bfld_value[2] << std::endl;
   AdvanceHelix(yInput, Bfld_value, hstep, yOut);

  // We are assuming a constant field: helix is exact
  //
  for( unsigned int i=0;i<nvar;i++)
  {
    yErr[i] = 0.0 ;
  }

  fBfieldValue=Bfld_value;
}

template <class Backend>
void
TemplateGUExactHelixStepper<Backend>
  ::StepWithoutErrorEstimate( const typename Backend::precision_v  yIn[],
                                    vecgeom::Vector3D<typename Backend::precision_v>   Bfld,
                                    typename Backend::precision_v  h,
                                    typename Backend::precision_v  yOut[]                  )
{
  // Assuming a constant field: solution is a helix

  AdvanceHelix(yIn, Bfld, h, yOut);

  std::cerr<<"TemplateGUExactHelixStepper::StepWithoutErrorEstimate"
           << "should *NEVER* be called. StepWithErrorEstimate must do the work." << std::endl;
}  


// ---------------------------------------------------------------------------
template <class Backend>
typename Backend::precision_v
TemplateGUExactHelixStepper<Backend>
  ::DistChord() const 
{
  // Implementation : must check whether h/R >  pi  !!
  //   If( h/R <  pi)   DistChord=h/2*std::tan(Ang_curve/4)                <
  //   Else             DistChord=2*R_helix    -- approximate.  True value ~ diameter

  typedef typename Backend::precision_v Double_v;
  Double_v distChord;
  Double_v Ang_curve= this->GetAngCurve();

  Double_v multiplicationFactor = (1+vecgeom::VECGEOM_IMPL_NAMESPACE::cos(0.5*(twopi-Ang_curve)));

  vecgeom::MaskedAssign( Ang_curve <= pi, (1-vecgeom::VECGEOM_IMPL_NAMESPACE::cos(0.5*Ang_curve)), &multiplicationFactor);
  vecgeom::MaskedAssign( Ang_curve >= twopi, 2., &multiplicationFactor);


  // if (Ang_curve<=pi)
  // {
  //   distChord=GetRadHelix()*(1-std::cos(0.5*Ang_curve));
  // }
  // else if(Ang_curve<twopi)
  // {
  //   distChord=GetRadHelix()*(1+std::cos(0.5*(twopi-Ang_curve)));
  // }
  // else
  // {
  //   distChord=2.*GetRadHelix();  
  // }

  distChord = this->GetRadHelix()*multiplicationFactor;

  return distChord;
}   

#endif  /* TemplateGUExactHelixStepper_h */
