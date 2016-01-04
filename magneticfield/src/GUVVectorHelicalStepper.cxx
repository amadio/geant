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
GUVVectorHelicalStepper::AdvanceHelix( const Vc::Vector<double>  yIn[],
                                             vecgeom::Vector3D<Vc::Vector<double> > Bfld,    
                                             double  h,
                                             Vc::Vector<double>  yHelix[],
                                             Vc::Vector<double>  yHelix2[] )
{
  // const G4int    nvar = 6;
 
  // OLD  const double approc_limit = 0.05;
  // OLD  approc_limit = 0.05 gives max.error=x^5/5!=(0.05)^5/5!=2.6*e-9
  // NEW  approc_limit = 0.005 gives max.error=x^5/5!=2.6*e-14

  const double approc_limit = 0.005;


  typedef Vc::Vector<double> Double_v;
  typedef vecgeom::Vector3D<Vc::Vector<double> > ThreeVectorSimd;
  typedef typename vecgeom::kVc::bool_v Bool_v;

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
  Double_v        B_v_P = Vc::sqrt( 1 - B_d_P * B_d_P); // Fraction of P perp to B
  Double_v        Theta = R_1 * h; // * B_v_P;

  //else inside else
  Double_v Theta2 = Theta  * Theta;
  Double_v Theta3 = Theta2 * Theta;
  Double_v Theta4 = Theta2 * Theta2;
  SinT     = Theta - 1.0/6.0 * Theta3;
  CosT     = 1 - 0.5 * Theta2 + 1.0/24.0 * Theta4;

  //if inside else
  Bool_v trigCond = Vc::abs(Theta) > approc_limit;
  // vecCore::MaskedAssign(&SinT, trigCond, Vc::sin(Theta));
  vecCore__MaskedAssignFunc(&SinT, trigCond, Vc::sin(Theta));
  // vecCore::MaskedAssign(&CosT, trigCond, Vc::cos(Theta));
  vecCore__MaskedAssignFunc(&CosT, trigCond, Vc::cos(Theta));
  
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

  R_Helix =Vc::abs( ptan/(fUnitConstant  * fParticleCharge*Bmag));
     

  // for too small magnetic fields there is no curvature
  // (include momentum here) FIXME

  Bool_v noCurvatureCond = (Vc::abs(R_1) < 1e-10) || (Bmag<1e-12);
  veccore__MaskedAssignFunc( &Theta,   noCurvatureCond, 1. );
  veccore__MaskedAssignFunc( &R,       noCurvatureCond, h  );
  veccore__MaskedAssignFunc( $R_Helix, noCurvatureCond, 0. );

  SetAngCurve(Vc::abs(Theta));
  SetCurve(Vc::abs(R));
  SetRadHelix(R_Helix);
}


//
//  Use the midpoint method to get an error estimate and correction
//  modified from G4ClassicalRK4: W.Wander <wwc@mit.edu> 12/09/97
//

void
GUVVectorHelicalStepper::StepWithErrorEstimate( const Vc::Vector<double> yInput[],
                                                const Vc::Vector<double>*,
                                                      double hstep,
                                                      Vc::Vector<double> yOut[],
                                                      Vc::Vector<double> yErr[]  )
{  
   const int nvar = 6;

   Vc::Vector<double>  yTemp[7], yIn[7] ;

   typedef vecgeom::Vector3D<Vc::Vector<double> > ThreeVectorSimd;

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


Vc::Vector<double> 
GUVVectorHelicalStepper::DistChord() const 
{
  // Check whether h/R >  pi  !!
  // Method DistLine is good only for <  pi

  Vc::Vector<double> Ang=GetAngCurve();
  Vc::Vector<double> returnValue;

  vecCore__MaskedAssignFunc( &returnValue, Ang<=pi,             GetRadHelix()*(1-Vc::cos(0.5*Ang)) );
  vecCore__MaskedAssignFunc( &returnValue, Ang>pi && Ang<twopi, GetRadHelix()*(1+Vc::cos(0.5*(twopi-Ang))) ); 
  vecCore__MaskedAssignFunc( &returnValue, Ang>= twopi,         2*GetRadHelix() ); 

  return returnValue;

}
