// class G4NystromRK4
//
// Class description:
//
// Integrate the equations of the motion of a particle in a magnetic field
// using 4th Runge-Kutta-Nystrom method with errors estimation 
// (ATL-SOFT-PUB-2009-01)
// Current form can be used only for 'pure' magnetic field.
// Notes: 1) field must be time-independent.
//        2) time is not integrated
// 
// History:
// - Created: I.Gavrilenko   15.05.2009   (as G4AtlasRK4)
// - Adaptations:  J. Apostolakis  May-Nov 2009
// -------------------------------------------------------------------

#ifndef GXNYSTROMRK4_HH
#define GXNYSTROMRK4_HH

#include "GPTypeDef.h"
#include "GXEquationOfMotion.h"

struct GXNystromRK4
{
  GXEquationOfMotion *m_fEq;

  G4double      m_lastField[3];
  G4double      m_fldPosition[3];
  G4double      m_iPoint   [3];
  G4double      m_mPoint   [3];
  G4double      m_fPoint   [3];

};

extern "C" {

FQUALIFIER
void GXNystromRK4_Constructor(GXNystromRK4 *This, 
			      GXEquationOfMotion* EqRhs);

FQUALIFIER
void GXNystromRK4_Stepper(GXNystromRK4 *This,
			  const G4double P[],const G4double dPdS[],
			  G4double Step,G4double Po[],G4double Err[]);

FQUALIFIER
G4double GXNystromRK4_DistChord(GXNystromRK4 *This);

FQUALIFIER
void  GXNystromRK4_RightHandSide(GXNystromRK4 *This, 
				 const G4double P[], 
				 G4double dPdS[]);

FQUALIFIER
void GXNystromRK4_GetField (GXNystromRK4 *This, const G4double P[4]);

FQUALIFIER
G4double GXNystromRK4_TruncatedError(GXNystromRK4 *This,
				     G4double hstep,
				     const G4double yarrout[],
				     const G4double yerr_vec[] );

}

#endif  // GPNYSTROMRK4
