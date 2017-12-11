//
// GUExactHelixStepper:
//
//  Displace using helix - for field at initial location
//
//  As the field is assumed constant, an error is not calculated.
// 
//  Author: J. Apostolakis, 28 Jan 2005
//     Implementation adapted from ExplicitEuler of W.Wander 
// -------------------------------------------------------------------

#include <cfloat>
#include <iostream>

#include "Constants.h"
using Constants::pi;
using Constants::twopi;

#include <cfloat>

#include "GUExactHelixStepper.h"
// #include "GUPhysicalConstants.h"
// #include "ThreeVector.h"
// #include "GULineSection.h"

GUExactHelixStepper::GUExactHelixStepper(VScalarEquationOfMotion* EqRhs) // TMagFieldEquation
   : VScalarHelicalStepper(EqRhs, 1),  // "Order" = 1 - not really applicable
    fBfieldValue(DBL_MAX), // , DBL_MAX, DBL_MAX),
    fPtrMagEqOfMot(EqRhs)
{
  ;
}

GUExactHelixStepper::~GUExactHelixStepper() {} 

void
GUExactHelixStepper::StepWithErrorEstimate( const double yInput[],
                                            const double*, // dydx
                                                  double charge, 
                                                  double hstep,
                                                  double yOut[],
                                                  double yErr[]  )
{  
   const unsigned int nvar = 6;

   ThreeVector Bfld_value;

   MagFieldEvaluate(yInput, Bfld_value);
   // std::cout << " Exact Helix: B-field:  Bx = " << Bfld_value[0]
   //           << " By= " << Bfld_value[1] << " Bz= " << Bfld_value[2] << std::endl;
   AdvanceHelix(yInput, Bfld_value, charge, hstep, yOut);

  // We are assuming a constant field: helix is exact
  //
  for( unsigned int i=0;i<nvar;i++)
  {
    yErr[i] = 0.0 ;
  }

  fBfieldValue=Bfld_value;
}

void
GUExactHelixStepper::StepWithoutErrorEstimate( const double  yIn[],
                                        ThreeVector   Bfld,
                                               double charge, 
                                        double  h,
                                        double  yOut[])
{
  // Assuming a constant field: solution is a helix

  AdvanceHelix(yIn, Bfld, charge, h, yOut);

  std::cerr<<"GUExactHelixStepper::StepWithoutErrorEstimate"
           << "should *NEVER* be called. StepWithErrorEstimate must do the work." << std::endl;
}  


// ---------------------------------------------------------------------------

double
GUExactHelixStepper::DistChord( double /*charge*/ ) const 
{
  // Implementation : must check whether h/R >  pi  !!
  //   If( h/R <  pi)   DistChord=h/2*std::tan(Ang_curve/4)                <
  //   Else             DistChord=2*R_helix    -- approximate.  True value ~ diameter

  double distChord;
  double Ang_curve=GetAngCurve();

  if (Ang_curve<=pi)
  {
    distChord=GetRadHelix()*(1-std::cos(0.5*Ang_curve));
  }
  else if(Ang_curve<twopi)
  {
    distChord=GetRadHelix()*(1+std::cos(0.5*(twopi-Ang_curve)));
  }
  else
  {
    distChord=2.*GetRadHelix();  
  }

  return distChord;
}   
