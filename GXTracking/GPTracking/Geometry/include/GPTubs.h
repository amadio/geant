#ifndef GPTUBS_HH
#define GPTUBS_HH

#include "GPVSolid.h"
#include "GPThreeVectorList.h"
#include "GPVoxelLimits.h"
#include "GPAffineTransform.h"

struct GPTubs
{
  GPVSolid fSolid;
  G4double fRMin, fRMax, fDz, fSPhi, fDPhi;
  G4bool fPhiFullTube;

  G4double sinCPhi, cosCPhi, cosHDPhiOT, cosHDPhiIT,
    sinSPhi, cosSPhi, sinEPhi, cosEPhi;
};

extern "C" {

  FQUALIFIER
  void GPTubs_Constructor(GPTubs *This, G4double pRMin, G4double pRMax,
			  G4double pDz, G4double pSPhi, G4double pDPhi );


FQUALIFIER
G4bool GPTubs_CalculateExtent(GPTubs *This, 
			      const EAxis       pAxis,
			      GPVoxelLimits     pVoxelLimit,
			      GPAffineTransform pTransform,
			      G4double          *pMin, 
			      G4double          *pMax    );

FQUALIFIER
EInside GPTubs_Inside( GEOMETRYLOC GPTubs *This, 
		       GPThreeVector p );

FQUALIFIER
GPThreeVector GPTubs_SurfaceNormal( GEOMETRYLOC GPTubs *This,
				    GPThreeVector p);

FQUALIFIER
GPThreeVector GPTubs_ApproxSurfaceNormal( GEOMETRYLOC  GPTubs *This, 
					  GPThreeVector p );

FQUALIFIER
G4double GPTubs_DistanceToIn2( GEOMETRYLOC GPTubs *This,
			       GPThreeVector p,
			       GPThreeVector v  );

FQUALIFIER
G4double GPTubs_DistanceToIn( GEOMETRYLOC GPTubs *This,
			      GPThreeVector p );

FQUALIFIER
G4double GPTubs_DistanceToOut2( GEOMETRYLOC GPTubs *This,
				GPThreeVector p,
				GPThreeVector v,
				const G4bool calcNorm,
				G4bool *validNorm,
				GPThreeVector *n    );

FQUALIFIER
G4double GPTubs_DistanceToOut( GEOMETRYLOC GPTubs *This,
			       GPThreeVector p );

FQUALIFIER
GPThreeVectorList
GPTubs_CreateRotatedVertices(GPTubs *This, 
			     GPAffineTransform pTransform ) ;

FQUALIFIER
G4double GPTubs_GetInnerRadius (GPTubs *This);

FQUALIFIER
G4double GPTubs_GetOuterRadius (GPTubs *This);

FQUALIFIER
G4double GPTubs_GetZHalfLength (GPTubs *This);

FQUALIFIER
G4double GPTubs_GetStartPhiAngle (GPTubs *This);

FQUALIFIER
G4double GPTubs_GetDeltaPhiAngle (GPTubs *This);

FQUALIFIER 
void GPTubs_Initialize(GPTubs *This);

FQUALIFIER 
void GPTubs_InitializeTrigonometry(GPTubs *This);

FQUALIFIER 
void GPTubs_CheckSPhiAngle(GPTubs *This,G4double sPhi);

FQUALIFIER 
void GPTubs_CheckDPhiAngle(GPTubs *This,G4double dPhi);

FQUALIFIER 
void GPTubs_CheckPhiAngles(GPTubs *This, 
			   G4double sPhi, G4double dPhi);

FQUALIFIER
void GPTubs_SetInnerRadius (GPTubs *This,
			    G4double newRMin);

FQUALIFIER
void GPTubs_SetOuterRadius (GPTubs *This,
			    G4double newRMax);

FQUALIFIER
void GPTubs_SetZHalfLength (GPTubs *This,
			    G4double newDz);

FQUALIFIER
void GPTubs_SetStartPhiAngle (GPTubs *This,
			      G4double newSPhi, G4bool compute);

FQUALIFIER
void GPTubs_SetDeltaPhiAngle (G4double newDPhi);

FQUALIFIER
G4double GPTubs_GetRMin (GPTubs *This);

FQUALIFIER
G4double GPTubs_GetRMax (GPTubs *This);

FQUALIFIER
G4double GPTubs_GetDz (GPTubs *This);

FQUALIFIER
G4double GPTubs_GetSPhi (GPTubs *This);

FQUALIFIER
G4double GPTubs_GetDPhi (GPTubs *This);

}

#endif
