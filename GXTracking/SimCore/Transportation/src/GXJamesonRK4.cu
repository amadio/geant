#include "GXJamesonRK4.h"
#include "GPUtils.h"
#include "GPLineSection.h"

#include "stdio.h"

FQUALIFIER
void GXJamesonRK4_Constructor(GXJamesonRK4 *This,
			      GXEquationOfMotion* Equation)
{
  This->fEquation_Rhs = Equation;
}

FQUALIFIER
void GXJamesonRK4_DumbStepper( GXJamesonRK4 *This,
			       const G4double  yIn[],
			       const G4double  dydx[],
			       G4double  h,
			       G4double  yOut[])
{
  const G4int nvar = 6;

  G4int i;
  G4double  hh = h*0.25;

  // 4-stage Jameson-Schmidt-Turkel Scheme
  // AiAA 1981 - p. 1259 Numerical solutions of the Euler equations by finite 
  //                     volume methods using Runge-Kutta time-stepping schemes

  for(i=0;i<nvar;++i) yOut[i] = yIn[i] + hh*dydx[i] ;     

  G4double dydxt[nvar];
  for(int k = 3 ; k >=1 ; --k) {
    GXJamesonRK4_RightHandSide(This, yOut, dydxt) ; 
    hh = h/k;
    for(i=0;i<nvar;++i) {
      yOut[i] = yIn[i] + hh*dydxt[i] ;     
    }
  }
}  

// G4MagErrorStepper::Stepper
FQUALIFIER
void GXJamesonRK4_Stepper( GXJamesonRK4 *This, 
			   const G4double yInput[],
			   const G4double dydx[],
			   G4double hstep,
			   G4double yOutput[],
			   G4double yError []      )
{  
  const G4int nvar = 6;

  G4int i;
  // correction for Richardson Extrapolation for the step double RK4
  G4double correction = 1./15. ;
  
  //  Saving yInput because yInput and yOutput can be aliases for same array
  
  for(i=0;i<nvar;i++) This->yInitial[i]=yInput[i];
  
  G4double halfStep = hstep * 0.5; 
  
  // Do two half steps
  
  GXJamesonRK4_DumbStepper(This, This->yInitial, dydx,
                             halfStep, This->yMiddle);
  GXJamesonRK4_RightHandSide(This, This->yMiddle, This->dydxMid);    
  GXJamesonRK4_DumbStepper(This, This->yMiddle, This->dydxMid, 
                             halfStep, yOutput); 
  
  // Store midpoint, chord calculation
  for(i=0;i<3;i++) This->fMidPoint[i] = This->yMiddle[i];  
  
  // Do a full Step
  GXJamesonRK4_DumbStepper(This,This->yInitial, dydx, hstep, This->yOneStep);
  for(i=0;i<nvar;i++) {
    yError [i] = yOutput[i] - This->yOneStep[i] ;
    yOutput[i] += yError[i]*correction ;  // Provides accuracy increased
    // by 1 order via the Richardson Extrapolation  
  }
  
  for(i=0;i<3;i++) {
    This->fInitialPoint[i] = This->yInitial[i];  
    This->fFinalPoint[i] = yOutput[i];  
  }
}

FQUALIFIER
void GXJamesonRK4_RightHandSide(GXJamesonRK4 *This,
                                  const  G4double y[],
                                  G4double dydx[] )
{
  GXEquationOfMotion_RightHandSide(This->fEquation_Rhs, y, dydx);
}

FQUALIFIER
G4double GXJamesonRK4_DistChord(GXJamesonRK4 *This)
{
  // Estimate the maximum distance from the curve to the chord
  //
  //  We estimate this using the distance of the midpoint to 
  //  chord (the line between 
  // 
  //  Method below is good only for angle deviations < 2 pi, 
  //   This restriction should not a problem for the Runge cutta methods, 
  //   which generally cannot integrate accurately for large angle deviations.
  G4double distChord; 

  if (This->fInitialPoint[0] != This->fFinalPoint[0] &&
      This->fInitialPoint[1] != This->fFinalPoint[1] &&
      This->fInitialPoint[2] != This->fFinalPoint[2] ) {

    distChord = GXJamesonRK4_DistLine( This );
  }else{
    GPThreeVector VecMF = 
      GPThreeVector_create (This->fMidPoint[0] - This->fInitialPoint[0],
			    This->fMidPoint[1] - This->fInitialPoint[1],
			    This->fMidPoint[2] - This->fInitialPoint[2]);
    distChord = GPThreeVector_mag(VecMF);
  }

  return distChord;
}


FQUALIFIER
G4double GXJamesonRK4_DistLine(GXJamesonRK4 *This)
{
  GPThreeVector OtherPnt  = GPThreeVector_create(This->fMidPoint[0],
						 This->fMidPoint[1],
						 This->fMidPoint[2]);
  GPThreeVector PntA = GPThreeVector_create(This->fInitialPoint[0],
						 This->fInitialPoint[1],
						 This->fInitialPoint[2]);
  GPThreeVector PntB = GPThreeVector_create(This->fFinalPoint[0],
						 This->fFinalPoint[1],
						 This->fFinalPoint[2]);

  GPThreeVector VecAtoB = GPThreeVector_sub(PntB,PntA);
  G4double ABdistanceSq = GPThreeVector_mag2(VecAtoB) ;

  G4double       dist_sq = 9999.;  
  GPThreeVector  VecAZ;
  G4double sq_VecAZ, inner_prod, unit_projection ; 

  VecAZ= GPThreeVector_sub(OtherPnt,PntA);
  sq_VecAZ = GPThreeVector_mag2(VecAZ);

  inner_prod= GPThreeVector_dot(VecAtoB,VecAZ );
   
  if( ABdistanceSq != 0.0 )
  {
    unit_projection = inner_prod/ABdistanceSq;

    if( (0. <= unit_projection ) && (unit_projection <= 1.0 ) )
    {
      dist_sq= sq_VecAZ -  unit_projection * inner_prod ;
    }
    else
    {
      if( unit_projection < 0. ) // A is the closest point
      {
        dist_sq= sq_VecAZ;  
      }
      else                       // B is the closest point
      {
	GPThreeVector   EndpointB = GPThreeVector_add(PntA,
						      VecAtoB);
        GPThreeVector   VecBZ =     GPThreeVector_sub(OtherPnt,
						      EndpointB);
        dist_sq =  GPThreeVector_mag2(VecBZ);
      }
    }
  }
  else
  {
    dist_sq = GPThreeVector_mag2(GPThreeVector_sub(OtherPnt,
						   PntA)) ;   
  }  
  if( dist_sq < 0.0 ) dist_sq = 0.0 ;

  return sqrt(dist_sq) ;  
}

FQUALIFIER
G4double GXJamesonRK4_TruncatedError(GXJamesonRK4 *This,
				     G4double hstep,
				     const G4double yarrout[],
				     const G4double yerr_vec[] )
{
  G4double dyerr = 0.0;
  G4double dyerr_pos_sq, dyerr_mom_rel_sq;  
  G4double dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

  vel_mag_sq = (yarrout[3]*yarrout[3] + yarrout[4]*yarrout[4] +
  		 yarrout[5]*yarrout[5] );
  inv_vel_mag_sq = 1.0 / vel_mag_sq; 

  dyerr_pos_sq = ( yerr_vec[0]*yerr_vec[0]+yerr_vec[1]*yerr_vec[1]+
          yerr_vec[2]*yerr_vec[2]);
  dyerr_mom_sq = ( yerr_vec[3]*yerr_vec[3]+yerr_vec[4]*yerr_vec[4]+
          yerr_vec[5]*yerr_vec[5]);
  dyerr_mom_rel_sq = dyerr_mom_sq * inv_vel_mag_sq;

  if( dyerr_pos_sq > ( dyerr_mom_rel_sq * hstep*hstep ) )
  {
    dyerr = sqrt(dyerr_pos_sq);
  }
  else
  {
    dyerr = sqrt(dyerr_mom_rel_sq) * hstep;
  }

  return dyerr;
}

