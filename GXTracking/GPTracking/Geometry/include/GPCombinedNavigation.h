#ifndef GPNORMALNAVIGATION_HH
#define GPNORMALNAVIGATION_HH

//#include <iomanip>

#include "GPTypeDef.h"
#include "GPThreeVector.h"
#include "GPVPhysicalVolume.h"
#include "GPLogicalVolume.h"
#include "GPVSolid.h"
#include "GPNavigationHistory.h"
#include "GPVoxelNavigation.h"

struct GPCombinedNavigation
{
    G4bool fCheck; 
};

extern "C" {

FQUALIFIER 
G4double GPCombinedNavigation_ComputeStep(
		       GPVoxelNavigation *vox,
		       GPThreeVector localPoint,
		       GPThreeVector localDirection,
		       const G4double currentProposedStepLength,
		       G4double *newSafety,
		       GPNavigationHistory *history,
		       G4bool *validExitNormal,
		       GPThreeVector *exitNormal,
		       G4bool *exiting,
		       G4bool *entering,
		       GEOMETRYLOC GPVPhysicalVolume *(*pBlockedPhysical) );

FQUALIFIER
G4bool GPCombinedNavigation_LevelLocate(
			GPVoxelNavigation *vox,
			GPNavigationHistory* history,
			GEOMETRYLOC const GPVPhysicalVolume* blockedVol,
			GPThreeVector globalPoint,
			const GPThreeVector* globalDirection,
			const G4bool pLocatedOnEdge, 
			GPThreeVector *localPoint );
FQUALIFIER 
G4double GPCombinedNavigation_ComputeSafety(
			GPVoxelNavigation *vox,
			GPThreeVector localPoint,
			const GPNavigationHistory *history );

}
#endif
