#ifndef GPCons_HH
#define GPCons_HH

#include "GPVSolid.h"

struct GPCons
{
  GPVSolid fSolid;
  G4double fRmin1, fRmin2, fRmax1, fRmax2, fDz, fSPhi, fDPhi;

  G4double sinCPhi, cosCPhi, cosHDPhiOT, cosHDPhiIT,
           sinSPhi, cosSPhi, sinEPhi, cosEPhi;

  G4bool fPhiFullCone;

};

extern "C" {

FQUALIFIER
void GPCons_Constructor(GPCons *This, 
			G4double  pRmin1, G4double pRmax1,
			G4double  pRmin2, G4double pRmax2,
			G4double pDz,
			G4double pSPhi, G4double pDPhi);

FQUALIFIER
G4double GPCons_GetInnerRadiusMinusZ( const GPCons *This );

FQUALIFIER
G4double GPCons_GetOuterRadiusMinusZ( const GPCons *This );

FQUALIFIER
G4double GPCons_GetInnerRadiusPlusZ( const GPCons *This );

FQUALIFIER
G4double GPCons_GetOuterRadiusPlusZ( const GPCons *This );

FQUALIFIER
G4double GPCons_GetZHalfLength( const GPCons *This );

FQUALIFIER  
G4double GPCons_GetStartPhiAngle( const GPCons *This );

FQUALIFIER
G4double GPCons_GetDeltaPhiAngle( const GPCons *This );

FQUALIFIER 
void GPCons_InitializeTrigonometry( GPCons *This );

FQUALIFIER 
void GPCons_CheckSPhiAngle( GPCons *This, G4double sPhi );

FQUALIFIER 
void GPCons_CheckDPhiAngle( GPCons *This, G4double dPhi );

FQUALIFIER 
void GPCons_CheckPhiAngles( GPCons *This, G4double sPhi, G4double dPhi );

FQUALIFIER 
EInside GPCons_Inside( GEOMETRYLOC const GPCons *This, 
		       GPThreeVector p);

FQUALIFIER 
GPThreeVectorList 
GPCons_CreateRotatedVertices(const GPCons *This, 
			     const GPAffineTransform pTransform);
FQUALIFIER
G4bool GPCons_CalculateExtent(const GPCons *This,
			      const EAxis pAxis,
			      GPVoxelLimits pVoxelLimit,
			      GPAffineTransform pTransform,
			      G4double* pMin, G4double* pMax);

FQUALIFIER
GPThreeVector GPCons_ApproxSurfaceNormal( GEOMETRYLOC const GPCons *This, 
					  GPThreeVector p );

FQUALIFIER
GPThreeVector GPCons_SurfaceNormal( GEOMETRYLOC const GPCons *This, 
				    GPThreeVector p);

FQUALIFIER 
G4double GPCons_DistanceToIn2( GEOMETRYLOC const GPCons *This,
			       GPThreeVector p,
			       GPThreeVector v);

FQUALIFIER
G4double GPCons_DistanceToIn( GEOMETRYLOC const GPCons *This, 
			      GPThreeVector p);

FQUALIFIER 
G4double GPCons_DistanceToOut2( GEOMETRYLOC const GPCons *This,
				GPThreeVector p,
				GPThreeVector v,
				const G4bool calcNorm,
				G4bool *validNorm,
				GPThreeVector *n);

FQUALIFIER
G4double GPCons_DistanceToOut( GEOMETRYLOC const GPCons *This, 
			       GPThreeVector p);

}
#endif
