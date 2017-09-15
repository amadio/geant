//
// Created:  16.10.2015 J. Apostolakis
//  Based on G4MagHelicalStepper  - meant for testing other steppers
// --------------------------------------------------------------------

// #include "G4SystemOfUnits.hh"

#include "Units.h"     
using fieldUnits::meter;  //  Update to GeantV units ASAP
using fieldUnits::GeV;
using fieldUnits::tesla;

#include "Constants.h"
using Constants::pi;
using Constants::twopi;

// #include "G4PhysicalConstants.hh"

#include "GUVVectorHelicalStepper.h"

// #include "GULineSection.h"
#include "TVectorMagFieldEquation.h"

// given a purely magnetic field a better approach than adding a straight line
// (as in the normal runge-kutta-methods) is to add helix segments to the
// current position

// Constant for determining unit conversion when using normal as integrand.
//
const double GUVVectorHelicalStepper::fUnitConstant = 0.299792458*(GeV/(tesla*meter));

GUVVectorHelicalStepper::GUVVectorHelicalStepper(GUVVectorEquationOfMotion *EqRhs,
                                                 unsigned int order              )
   : GUVVectorIntegrationStepper(EqRhs, order, 6,  6), //integrate over 6 variables only, // state could be 8 - also t, E
                               
     fPtrMagEqOfMot(EqRhs), fAngCurve(0.), frCurve(0.), frHelix(0.), fParticleCharge(0.0)
{
}

GUVVectorHelicalStepper::~GUVVectorHelicalStepper()
{
}

void
GUVVectorHelicalStepper::AdvanceHelix( const Double_v  yIn[],
                                             ThreeVectorSimd Bfld,
                                       const Double_v& charge,
                                       const Double_v& h,
                                             Double_v  yHelix[],
                                             Double_v  yHelix2[] )
{
  // const G4int    nvar = 6;
 
  // OLD  const double approc_limit = 0.05;
  // OLD  approc_limit = 0.05 gives max.error=x^5/5!=(0.05)^5/5!=2.6*e-9
  // NEW  approc_limit = 0.005 gives max.error=x^5/5!=2.6*e-14

  const double approc_limit = 0.005;

  // using vecCore::MaskedAssign;

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
  LinearStep( yIn, h, yHelix );
  
  //from else statement
  ThreeVectorSimd Bnorm = (1.0/Bmag)*Bfld;
  ThreeVectorSimd B_x_P = Bnorm.Cross(initTangent);
  Double_v        B_d_P = Bnorm.Dot(initTangent); // this is the fraction of P parallel to B
  ThreeVectorSimd vpar  = B_d_P * Bnorm;       // the component parallel      to B
  ThreeVectorSimd vperp = initTangent - vpar;  // the component perpendicular to B
  Double_v        B_v_P = vecCore::math::Sqrt( 1 - B_d_P * B_d_P); // Fraction of P perp to B
  Double_v        Theta = R_1 * h; // * B_v_P;

  //else inside else
  Double_v Theta2 = Theta  * Theta;
  Double_v Theta3 = Theta2 * Theta;
  Double_v Theta4 = Theta2 * Theta2;
  SinT     = Theta - 1.0/6.0 * Theta3;
  CosT     = 1 - 0.5 * Theta2 + 1.0/24.0 * Theta4;

  //if inside else
  Bool_v trigCond = vecCore::math::Abs(Theta) > approc_limit;
  // vecCore::MaskedAssign( SinT, trigCond, vecCore::math::Sin(Theta));
  vecCore__MaskedAssignFunc(SinT, trigCond, vecCore::math::Sin(Theta));
  // vecCore::MaskedAssign( CosT, trigCond, vecCore::math::Cos(Theta));
  vecCore__MaskedAssignFunc(CosT, trigCond, vecCore::math::Cos(Theta));
  
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

  using vecCore::math::Abs;

  R_Helix =Abs( ptan/(fUnitConstant  * charge*Bmag));
  // Was: 
  //    R_Helix =vecCore::math::Abs( ptan/(fUnitConstant  * charge*Bmag));
  
  // for too small magnetic fields there is no curvature
  // (include momentum here) FIXME
  // Bool_v
  vecCore::Mask<Double_v>
     noCurvatureCond = ( Abs(R_1) < 1e-10) || (Bmag<1e-12);

  vecCore::MaskedAssign( Theta,   noCurvatureCond, Double_v(1.) );
  // vecCore__MaskedAssignFunc( Theta,   noCurvatureCond, 1. );  
  // vecCore::MaskedAssign( R,       noCurvatureCond, h  );
  vecCore::MaskedAssign( R,       noCurvatureCond, Double_v(h)  );  
  // vecCore::MaskedAssign( R_Helix, noCurvatureCond, 0. );
  vecCore::MaskedAssign( R_Helix, noCurvatureCond, Double_v(0.) );

  // R       = vecCore::Blend( noCurvatureCond, h  );  
  // R_Helix = vecCore::Blend( noCurvatureCond, 0, Abs( ptan/(fUnitConstant  * charge*Bmag)) );
  
  SetAngCurve( vecCore::math::Abs(Theta) );
  SetCurve(    vecCore::math::Abs(R) );
  SetRadHelix( R_Helix    );
}


//
//  Use the midpoint method to get an error estimate and correction
//  modified from G4ClassicalRK4: W.Wander <wwc@mit.edu> 12/09/97
//

void
GUVVectorHelicalStepper::StepWithErrorEstimate( const Double_v yInput[],
                                                const Double_v*,         // dydx
                                                const Double_v& charge,
                                                const Double_v& hstep,
                                                      Double_v yOut[],
                                                      Double_v yErr[]  )
{
   const int Nvar = 6;

   Double_v  yMid[Nvar], yIn[Nvar] ;
   Double_v  ySingleStep[Nvar];

   ThreeVectorSimd Bfld_initial, Bfld_midpoint;
   
   //  Saving yInput because yInput and yOut can be aliases for same array

   for(unsigned int i=0;i<Nvar;i++) { yIn[i]=yInput[i]; }

   Double_v halfStep = hstep * 0.5;

   MagFieldEvaluate(yIn, Bfld_initial) ;

   // Do first half step

   StepWithoutErrorEstimate(yIn, Bfld_initial,  charge, halfStep, yMid);
   MagFieldEvaluate(yMid, Bfld_midpoint) ;

   // Do a full Step
   StepWithoutErrorEstimate(yIn, Bfld_initial,  charge, hstep,    ySingleStep);
   
   // Do second half step   
   StepWithoutErrorEstimate(yMid, Bfld_midpoint, charge, halfStep, yOut); 

   // Error estimation
   for(unsigned int i=0; i<Nvar; ++i)
   {
     yErr[i] = yOut[i] - ySingleStep[i] ;
   }
   
   return;
}


Geant::Double_v 
GUVVectorHelicalStepper::DistChord() const 
{
  // Check whether h/R >  pi  !!
  // Method DistLine is good only for <  pi

  Double_v Ang=GetAngCurve();
  Double_v returnValue;

  vecCore__MaskedAssignFunc( returnValue, Ang<=pi,             GetRadHelix()*(1-vecCore::math::Cos(0.5*Ang)) );
  vecCore__MaskedAssignFunc( returnValue, Ang>pi && Ang<twopi, GetRadHelix()*(1+vecCore::math::Cos(0.5*(twopi-Ang))) ); 
  vecCore__MaskedAssignFunc( returnValue, Ang>= twopi,         2*GetRadHelix() ); 

  return returnValue;

}
