#ifndef GPTRD_HH
#define GPTRD_HH

#include "GPTypeDef.h"
#include "GPGeomdefs.h"
#include "GPThreeVector.h"
#include "GPThreeVectorList.h"

#include "GPVSolid.h"
#include "GPVoxelLimits.h"
#include "GPAffineTransform.h"

struct GPTrd
{
  GPVSolid fSolid;
  G4double fDx1,fDx2,fDy1,fDy2,fDz;
};

extern "C" {

FQUALIFIER
void GPTrd_Constructor(GPTrd *This,
		       G4double pdx1,  G4double pdx2,
		       G4double pdy1,  G4double pdy2,
		       G4double pdz );

FQUALIFIER
void GPTrd_CheckAndSetAllParameters ( GPTrd *This,
				      G4double pdx1,  G4double pdx2,
				      G4double pdy1,  G4double pdy2,
				      G4double pdz );

FQUALIFIER
void GPTrd_SetAllParameters (GPTrd *This, 
			     G4double pdx1, G4double pdx2, G4double pdy1, 
			     G4double pdy2, G4double pdz );

FQUALIFIER
G4bool GPTrd_CalculateExtent(GPTrd *This, 
			     const EAxis pAxis,
			     GPVoxelLimits pVoxelLimit,
			     GPAffineTransform pTransform,
			     G4double *pMin, G4double *pMax );

FQUALIFIER
EInside GPTrd_Inside( GEOMETRYLOC GPTrd *This,
		      GPThreeVector p );

FQUALIFIER
GPThreeVector GPTrd_SurfaceNormal( GEOMETRYLOC GPTrd *This,
				   GPThreeVector p );

FQUALIFIER
GPThreeVector GPTrd_ApproxSurfaceNormal( GEOMETRYLOC GPTrd *This,
					 GPThreeVector p );

FQUALIFIER
G4double GPTrd_DistanceToIn2( GEOMETRYLOC GPTrd *This, 
			      GPThreeVector p,
			      GPThreeVector v );

FQUALIFIER
G4double GPTrd_DistanceToIn( GEOMETRYLOC GPTrd *This,
			     GPThreeVector p );

FQUALIFIER
G4double GPTrd_DistanceToOut2( GEOMETRYLOC GPTrd *This,
			       GPThreeVector p,
			       GPThreeVector v,
			       const G4bool calcNorm,
			       G4bool *validNorm,
			       GPThreeVector *n );

FQUALIFIER
G4double GPTrd_DistanceToOut( GEOMETRYLOC GPTrd *This, 
			      GPThreeVector p );

FQUALIFIER
GPThreeVectorList
GPTrd_CreateRotatedVertices(GPTrd *This, 
			    GPAffineTransform pTransform );

FQUALIFIER
G4double GPTrd_GetXHalfLength1(GPTrd *This);

FQUALIFIER
G4double GPTrd_GetXHalfLength2(GPTrd *This);

FQUALIFIER
G4double GPTrd_GetYHalfLength1(GPTrd *This);

FQUALIFIER
G4double GPTrd_GetYHalfLength2(GPTrd *This);

FQUALIFIER
G4double GPTrd_GetZHalfLength(GPTrd *This);

FQUALIFIER
void GPTrd_SetXHalfLength1(GPTrd *This,
			   G4double val);

FQUALIFIER
void GPTrd_SetXHalfLength2(GPTrd *This,
			   G4double val);

FQUALIFIER
void GPTrd_SetYHalfLength1(GPTrd *This,
			   G4double val);

FQUALIFIER
void GPTrd_SetYHalfLength2(GPTrd *This,
			   G4double val);

FQUALIFIER
void GPTrd_SetZHalfLength(GPTrd *This,
			  G4double val);

}

#endif
