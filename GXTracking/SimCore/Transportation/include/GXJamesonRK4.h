#ifndef GXJamesonRK4_HH
#define GXJamesonRK4_HH

#include "GPTypeDef.h"
#include "GPThreeVector.h"
#include "GXEquationOfMotion.h"
#include "GXJamesonRK4.h"

struct GXJamesonRK4
{
  GXEquationOfMotion *fEquation_Rhs;
  // STATE
  G4double yInitial[6], yMiddle[6], dydxMid[6], yOneStep[6];
  G4double fInitialPoint[3], fMidPoint[3], fFinalPoint[3];

};

extern "C" {

FQUALIFIER
void GXJamesonRK4_Constructor( GXJamesonRK4 *This, 
			       GXEquationOfMotion* EqRhs);

FQUALIFIER
void GXJamesonRK4_DumbStepper( GXJamesonRK4 *This,
			       const G4double  yIn[],
			       const G4double  dydx[],
			       G4double  h,
			       G4double  yOut[]);

FQUALIFIER
void GXJamesonRK4_Stepper( GXJamesonRK4 *This, 
			   const G4double yInput[],
			   const G4double dydx[],
			   G4double hstep,
			   G4double yOutput[],
			   G4double yError [] );

FQUALIFIER 
void GXJamesonRK4_RightHandSide(GXJamesonRK4 *This,
				const  G4double y[], 
				G4double dydx[] );   

FQUALIFIER
G4double GXJamesonRK4_DistChord(GXJamesonRK4 *This);

FQUALIFIER
G4double GXJamesonRK4_DistLine(GXJamesonRK4 *This);

FQUALIFIER
G4double GXJamesonRK4_TruncatedError(GXJamesonRK4 *This,
				     G4double hstep,
				     const G4double yarrout[],
				     const G4double yerr_vec[]);
}

#endif
