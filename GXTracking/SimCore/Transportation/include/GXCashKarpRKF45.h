#ifndef GXCashKarpRKF45_HH
#define GXCashKarpRKF45_HH

#include "GPTypeDef.h"
#include "GPLineSection.h"

#include "GXEquationOfMotion.h"

struct GXCashKarpRKF45
{
  GXEquationOfMotion *fEquation_Rhs;

  G4double fLastStepLength;
  G4double fLastInitialVector[6], fLastFinalVector[6],
    fLastDyDx[6], fMidVector[6], fMidError[6];
};

extern "C" {

FQUALIFIER
void GXCashKarpRKF45_Constructor( GXCashKarpRKF45 *This, 
				  GXEquationOfMotion* EqRhs);

FQUALIFIER
void GXCashKarpRKF45_Stepper( GXCashKarpRKF45 *This, 
			      const G4double yInput[],
			      const G4double dydx[],
			      G4double Step,
			      G4double yOut[],
			      G4double yErr[]);

FQUALIFIER
G4double GXCashKarpRKF45_DistChord(GXCashKarpRKF45 *This);

FQUALIFIER 
void GXCashKarpRKF45_RightHandSide(GXCashKarpRKF45 *This,
				   const  G4double y[], 
				   G4double dydx[] );

FQUALIFIER
G4double GXCashKarpRKF45_TruncatedError(GXCashKarpRKF45 *This,
					G4double hstep,
					const G4double yarrout[],
					const G4double yerr_vec[] );
}

#endif
