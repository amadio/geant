#ifndef GXClassicalRK4_HH
#define GXClassicalRK4_HH

#include "GPTypeDef.h"
#include "GPThreeVector.h"
#include "GXEquationOfMotion.h"
#include "GXClassicalRK4.h"

struct GXClassicalRK4
{
  GXEquationOfMotion *fEquation_Rhs;
  // STATE
  G4double yInitial[6], yMiddle[6], dydxMid[6], yOneStep[6];
  G4double fInitialPoint[3], fMidPoint[3], fFinalPoint[3];

};

extern "C" {

FQUALIFIER
void GXClassicalRK4_Constructor( GXClassicalRK4 *This, 
				 GXEquationOfMotion* EqRhs);

FQUALIFIER
void GXClassicalRK4_DumbStepper( GXClassicalRK4 *This,
				 const G4double  yIn[],
				 const G4double  dydx[],
				 G4double  h,
				 G4double  yOut[]);

FQUALIFIER
void GXClassicalRK4_Stepper( GXClassicalRK4 *This, 
			     const G4double yInput[],
			     const G4double dydx[],
			     G4double hstep,
			     G4double yOutput[],
			     G4double yError [] );

FQUALIFIER 
void GXClassicalRK4_RightHandSide(GXClassicalRK4 *This,
				  const  G4double y[], 
				  G4double dydx[] );   

FQUALIFIER
G4double GXClassicalRK4_DistChord(GXClassicalRK4 *This);

FQUALIFIER
G4double GXClassicalRK4_DistLine(GXClassicalRK4 *This);

FQUALIFIER
G4double GXClassicalRK4_TruncatedError(GXClassicalRK4 *This,
                                       G4double hstep,
				       const G4double yarrout[],
				       const G4double yerr_vec[]);
}

#endif
