#ifndef GPAuxiliaryNavServices_HH
#define GPAuxiliaryNavServices_HH

#include "GPTypeDef.h"
#include "GPThreeVector.h"
#include "GPAffineTransform.h"
#include "GPVSolid.h"

extern "C" {

FQUALIFIER 
G4bool 
GPAuxiliaryNavServices_CheckPointOnSurface( GEOMETRYLOC const GPVSolid* sampleSolid, 
					    GPThreeVector localPoint, 
					    GPThreeVector *globalDirection, 
					    GPAffineTransform sampleTransform,
					    const G4bool locatedOnEdge);
  
FQUALIFIER 
G4bool 
GPAuxiliaryNavServices_CheckPointExiting( GEOMETRYLOC const GPVSolid* sampleSolid, 
					  const GPThreeVector localPoint, 
					  const GPThreeVector* globalDirection, 
					  const GPAffineTransform sampleTransform );

}

#endif
