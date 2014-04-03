#ifndef GPOrb_HH
#define GPOrb_HH

#include "GPVSolid.h"
#include "GPGeomdefs.h"

struct GPOrb
{
  GPVSolid fSolid;
  G4double fRmax;
  G4double fRmaxTolerance;
};

extern "C" {

  FQUALIFIER
  void GPOrb_Constructor(GPOrb *This, G4double radius);

FQUALIFIER
G4bool GPOrb_CalculateExtent(const GPOrb *This,
			     const EAxis pAxis,
			     GPVoxelLimits pVoxelLimit,
			     GPAffineTransform pTransform,
			     G4double* pMin, G4double* pMax);

FQUALIFIER
EInside GPOrb_Inside( GEOMETRYLOC const GPOrb *This, 
		      GPThreeVector p);

FQUALIFIER
GPThreeVector GPOrb_SurfaceNormal( GEOMETRYLOC const GPOrb *This, 
				   GPThreeVector p);

FQUALIFIER
G4double GPOrb_DistanceToIn2( GEOMETRYLOC const GPOrb *This, 
			      GPThreeVector p,
			      GPThreeVector v);

FQUALIFIER
G4double GPOrb_DistanceToIn( GEOMETRYLOC const GPOrb *This, 
			     GPThreeVector p);

FQUALIFIER
G4double GPOrb_DistanceToOut2( GEOMETRYLOC const GPOrb *This, 
			       GPThreeVector p,GPThreeVector v,
                               const G4bool calcNorm,
			       G4bool *validNorm,GPThreeVector *n);

FQUALIFIER
G4double GPOrb_DistanceToOut( GEOMETRYLOC const GPOrb *This, 
			      GPThreeVector p );

}

#endif
