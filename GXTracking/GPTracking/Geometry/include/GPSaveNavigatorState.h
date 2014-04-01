#ifndef GPSaveNavigatorState_HH
#define GPSaveNavigatorState_HH 1

#include "GPThreeVector.h"
#include "GPVPhysicalVolume.h"

struct GPSaveNavigatorState
{ 
  GPThreeVector sExitNormal;  
  G4bool sValidExitNormal;    
  G4bool sEntering, sExiting;
  GPVPhysicalVolume* spBlockedPhysicalVolume;
  G4int sBlockedReplicaNo;  
  G4int sLastStepWasZero; 
  
  //  Potentially relevant
  //
  G4bool sLocatedOutsideWorld;
  GPThreeVector sLastLocatedPointLocal; 
  G4bool sEnteredDaughter, sExitedMother;
  GPThreeVector  sPreviousSftOrigin;
  G4double       sPreviousSafety; 
};

#endif
