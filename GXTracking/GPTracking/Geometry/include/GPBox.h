#ifndef GPBOX_HH
#define GPBOX_HH

#include "GPVSolid.h"
#include "GPThreeVectorList.h"
#include "GPVoxelLimits.h"
#include "GPAffineTransform.h"

struct GPBox
{
  GPVSolid fSolid;
  G4double fDx,fDy,fDz;
};

extern "C" {

FQUALIFIER
  void GPBox_Constructor(GPBox *This, G4double x, G4double y, G4double z);

FQUALIFIER
void GPBox_SetXHalfLength(GPBox *This,
                          G4double dx);

FQUALIFIER
void GPBox_SetYHalfLength(GPBox *This,
                          G4double dy);

FQUALIFIER
void GPBox_SetZHalfLength(GPBox *This,
                          G4double dz);

FQUALIFIER
G4bool GPBox_CalculateExtent(GPBox *This,
                             const EAxis pAxis,
                             GPVoxelLimits pVoxelLimit,
                             GPAffineTransform pTransform,
                             G4double *pMin, G4double *pMax);

FQUALIFIER
EInside GPBox_Inside( GEOMETRYLOC GPBox *This,
		      GPThreeVector p);

FQUALIFIER
GPThreeVector GPBox_SurfaceNormal( GEOMETRYLOC GPBox *This,
				   GPThreeVector p);

FQUALIFIER
GPThreeVector GPBox_ApproxSurfaceNormal( GEOMETRYLOC GPBox *This,
					 GPThreeVector p );

FQUALIFIER
G4double GPBox_DistanceToIn2( GEOMETRYLOC GPBox *This,
			      GPThreeVector p,
			      GPThreeVector v);

FQUALIFIER
G4double GPBox_DistanceToIn( GEOMETRYLOC GPBox *This, 
			     GPThreeVector p);
  
  FQUALIFIER
G4double GPBox_DistanceToOut2( GEOMETRYLOC GPBox *This,
			       GPThreeVector p,
			       GPThreeVector v,
			       const G4bool calcNorm,
			       G4bool *validNorm, GPThreeVector *n);

FQUALIFIER
G4double GPBox_DistanceToOut( GEOMETRYLOC GPBox *This,
			      GPThreeVector p);

FQUALIFIER
GPThreeVectorList
GPBox_CreateRotatedVertices(GPBox *This,
                            GPAffineTransform pTransform);

FQUALIFIER
G4double GPBox_GetXHalfLength(GPBox *This);

FQUALIFIER    
G4double GPBox_GetYHalfLength(GPBox *This);

FQUALIFIER
G4double GPBox_GetZHalfLength(GPBox *This);

}

#endif
