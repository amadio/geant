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

#include "GUVHelicalStepper.h"

// #include "GULineSection.h"
#include "TMagFieldEquation.h"

// given a purely magnetic field a better approach than adding a straight line
// (as in the normal runge-kutta-methods) is to add helix segments to the
// current position

// Constant for determining unit conversion when using normal as integrand.
//
const double GUVHelicalStepper::fUnitConstant = 0.299792458*(GeV/(tesla*meter));

GUVHelicalStepper::GUVHelicalStepper(GUVEquationOfMotion *EqRhs,
                                     unsigned int order)
   : GUVIntegrationStepper(EqRhs,
                           order,
                           6,   // integrate over 6 variables only !! ( position & momentum )
                           6),  // state could be 8 - also t, E
                               
     fPtrMagEqOfMot(EqRhs), fAngCurve(0.), frCurve(0.), frHelix(0.)
{
}

GUVHelicalStepper::~GUVHelicalStepper()
{
}

void
GUVHelicalStepper::AdvanceHelix( const double  yIn[],
                                 ThreeVector   Bfld,
                                       double  charge, 
                                       double  h,
                                       double  yHelix[],
                                       double  yHelix2[] )
{
  // const G4int    nvar = 6;
 
  // OLD  const double approc_limit = 0.05;
  // OLD  approc_limit = 0.05 gives max.error=x^5/5!=(0.05)^5/5!=2.6*e-9
  // NEW  approc_limit = 0.005 gives max.error=x^5/5!=2.6*e-14

  const double approc_limit = 0.005;
  ThreeVector  Bnorm, B_x_P, vperp, vpar;

  double B_d_P;
  double B_v_P;
  double Theta;
  double R_1;
  double R_Helix;
  double CosT2, SinT2, CosT, SinT;
  ThreeVector positionMove, endTangent;

  double Bmag = Bfld.Mag2();
  const double *pIn = yIn+3;
  ThreeVector initVelocity= ThreeVector( pIn[0], pIn[1], pIn[2]);
  double      velocityVal = initVelocity.Mag2();
  ThreeVector initTangent = (1.0/velocityVal) * initVelocity;
  
  R_1=GetInverseCurve(velocityVal,charge,Bmag);

  // for too small magnetic fields there is no curvature
  // (include momentum here) FIXME

  if( (std::fabs(R_1) < 1e-10)||(Bmag<1e-12) )
  {
    LinearStep( yIn, h, yHelix );

    // Store and/or calculate parameters for chord distance

    SetAngCurve(1.);     
    SetCurve(h);
    SetRadHelix(0.);
  }
  else
  {
    Bnorm = (1.0/Bmag)*Bfld;

    // calculate the direction of the force
    
    B_x_P = Bnorm.Cross(initTangent);

    // parallel and perp vectors

    B_d_P = Bnorm.Dot(initTangent); // this is the fraction of P parallel to B

    vpar = B_d_P * Bnorm;       // the component parallel      to B
    vperp= initTangent - vpar;  // the component perpendicular to B
    
    B_v_P  = std::sqrt( 1 - B_d_P * B_d_P); // Fraction of P perp to B

    // calculate  the stepping angle
  
    Theta   = R_1 * h; // * B_v_P;

    // Trigonometrix
      
    if( std::fabs(Theta) > approc_limit )
    {
       SinT     = std::sin(Theta);
       CosT     = std::cos(Theta);
    }
    else
    {
      double Theta2 = Theta*Theta;
      double Theta3 = Theta2 * Theta;
      double Theta4 = Theta2 * Theta2;
      SinT     = Theta - 1.0/6.0 * Theta3;
      CosT     = 1 - 0.5 * Theta2 + 1.0/24.0 * Theta4;
    }

    // the actual "rotation"

    double R = 1.0 / R_1;

    positionMove  = R * ( SinT * vperp + (1-CosT) * B_x_P) + h * vpar;
    endTangent    = CosT * vperp + SinT * B_x_P + vpar;

    // Store the resulting position and tangent

    yHelix[0]   = yIn[0] + positionMove.x(); 
    yHelix[1]   = yIn[1] + positionMove.y(); 
    yHelix[2]   = yIn[2] + positionMove.z();
    yHelix[3] = velocityVal * endTangent.x();
    yHelix[4] = velocityVal * endTangent.y();
    yHelix[5] = velocityVal * endTangent.z();

    // Store 2*h step Helix if exist

    if(yHelix2)
    {
      SinT2     = 2.0 * SinT * CosT;
      CosT2     = 1.0 - 2.0 * SinT * SinT;
      endTangent    = (CosT2 * vperp + SinT2 * B_x_P + vpar);
      positionMove  = R * ( SinT2 * vperp + (1-CosT2) * B_x_P) + h*2 * vpar;
      
      yHelix2[0]   = yIn[0] + positionMove.x(); 
      yHelix2[1]   = yIn[1] + positionMove.y(); 
      yHelix2[2]   = yIn[2] + positionMove.z(); 
      yHelix2[3] = velocityVal * endTangent.x();
      yHelix2[4] = velocityVal * endTangent.y();
      yHelix2[5] = velocityVal * endTangent.z();
    }

    // Store and/or calculate parameters for chord distance

    double ptan=velocityVal*B_v_P;

    // double particleCharge = fPtrMagEqOfMot->FCof() / (eplus*c_light); 
    R_Helix =std::abs( ptan/(fUnitConstant  * charge * Bmag));
       
    SetAngCurve(std::abs(Theta));
    SetCurve(std::abs(R));
    SetRadHelix(R_Helix);
  }
}

//
//  Use the midpoint method to get an error estimate and correction
//  modified from G4ClassicalRK4: W.Wander <wwc@mit.edu> 12/09/97
//

void
GUVHelicalStepper::StepWithErrorEstimate( const double yInput[],
                                          const double*,      // dydx: Not relevant
                                          const double charge,
                                          double hstep,
                                          double yOut[],
                                          double yErr[]  )
{  
   const int nvar = 6;

   double      yTemp[7], yIn[7] ;
   ThreeVector Bfld_initial, Bfld_midpoint;
   
   //  Saving yInput because yInput and yOut can be aliases for same array

   for(unsigned int i=0;i<nvar;i++) { yIn[i]=yInput[i]; }

   double h = hstep * 0.5; 

   MagFieldEvaluate(yIn, Bfld_initial) ;      

   // Do two half steps

   StepWithoutErrorEstimate(yIn,   Bfld_initial,  charge, h, yTemp);
   MagFieldEvaluate(yTemp, Bfld_midpoint) ;     
   StepWithoutErrorEstimate(yTemp, Bfld_midpoint, charge, h, yOut); 

   // Do a full Step

   h = hstep ;
   StepWithoutErrorEstimate(yIn,   Bfld_initial,  charge, h, yTemp);

   // Error estimation

   for(unsigned int i=0; i<nvar; ++i)
   {
     yErr[i] = yOut[i] - yTemp[i] ;
   }
   
   return;
}

double
GUVHelicalStepper::DistChord(double /*charge*/ ) const 
{
  // Check whether h/R >  pi  !!
  // Method DistLine is good only for <  pi

  double Ang=GetAngCurve();
  if(Ang<=pi)
  {
    return GetRadHelix()*(1-std::cos(0.5*Ang));
  }
  else
  {
    if(Ang<twopi)
    {
      return GetRadHelix()*(1+std::cos(0.5*(twopi-Ang)));
    }
    else  // return Diameter of projected circle
    {
      return 2*GetRadHelix();
    }
  }
}
