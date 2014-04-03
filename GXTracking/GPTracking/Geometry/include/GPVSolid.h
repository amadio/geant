#ifndef GPVSOLID_HH
#define GPVSOLID_HH

#include "GPGeomdefs.h"
#include "GPThreeVector.h"
#include "GPThreeVectorList.h"
#include "GPVoxelLimits.h"
#include "GPAffineTransform.h"

//enum ESolid { kBox, kTrd };
enum ESolid { kBox, kOrb, kTubs, kCons, kPolyCone, kTrd };

struct GPVSolid
{
  ESolid fType;
};

extern "C" {

FQUALIFIER
void GPVSolid_Constructor(GPVSolid *This,
                          ESolid kType);

FQUALIFIER
void GPVSolid_ClipCrossSection( GPThreeVectorList* pVertices,
                                const G4int pSectionIndex,
                                GPVoxelLimits* pVoxelLimit,
                                const EAxis pAxis, 
				G4double *pMin, G4double *pMax);

FQUALIFIER
void GPVSolid_ClipBetweenSections( GPThreeVectorList* pVertices,
                                   const G4int pSectionIndex,
                                   GPVoxelLimits* pVoxelLimit,
                                   const EAxis pAxis, 
				   G4double *pMin, G4double *pMax);

FQUALIFIER
void
GPVSolid_CalculateClippedPolygonExtent(GPThreeVectorList pPolygon,
				       GPVoxelLimits* pVoxelLimit,
				       const EAxis pAxis, 
				       G4double *pMin,
				       G4double *pMax);

FQUALIFIER
void GPVSolid_ClipPolygon( GPThreeVectorList pPolygon,
                           GPVoxelLimits* pVoxelLimit,
                           const EAxis  ); 

FQUALIFIER
void
GPVSolid_ClipPolygonToSimpleLimits( GPThreeVectorList pPolygon,
				    GPThreeVectorList* outputPolygon,
				    GPVoxelLimits* pVoxelLimit       );

FQUALIFIER
EInside GPVSolid_Inside( GEOMETRYLOC const GPVSolid *This, 
			 GPThreeVector p);

FQUALIFIER
GPThreeVector GPVSolid_SurfaceNormal( GEOMETRYLOC const GPVSolid *This, 
				      GPThreeVector p);

FQUALIFIER
G4double GPVSolid_DistanceToIn( GEOMETRYLOC const GPVSolid *This, 
				GPThreeVector p);
  
FQUALIFIER
G4double GPVSolid_DistanceToOut( GEOMETRYLOC const GPVSolid *This, 
				 GPThreeVector p);
  
FQUALIFIER
G4double GPVSolid_DistanceToIn2( GEOMETRYLOC const GPVSolid *This, 
				 GPThreeVector p,
				 GPThreeVector v);

FQUALIFIER
G4double GPVSolid_DistanceToOut2( GEOMETRYLOC const GPVSolid *This,
				  GPThreeVector p,
				  GPThreeVector v,
				  const G4bool calcNorm,
				  G4bool *validNorm,
				  GPThreeVector *n);

FQUALIFIER
G4bool GPVSolid_CalculateExtent( const GPVSolid *This,
				 const EAxis pAxis,
				 GPVoxelLimits pVoxelLimit,
				 GPAffineTransform pTransform,
				 G4double* pMin, G4double* pMax);

}

#endif
