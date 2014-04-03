#ifndef GXeTrackingManager_H
#define GXeTrackingManager_H 1

#include "GPTypeDef.h"
#include "GXeProcessManager.h"       
#include "GXeSteppingManager.h"      
#include "GXTrack.h"                

class GXeTrackingManager 
{
public:

  FQUALIFIER GXeTrackingManager(GXeProcessManager* aProcessManager,
				GXeSteppingManager* aSteppingManager);

  FQUALIFIER ~GXeTrackingManager() {};

  FQUALIFIER GXTrack* GetTrack();
  FQUALIFIER GXeSteppingManager* GetSteppingManager();

  FQUALIFIER void SetNumberOfSteps(G4int aValue);
  FQUALIFIER void ProcessGPIL(GXTrack* apValue);
  FQUALIFIER void ProcessDoIt(GXTrack* apValue);

private:

  G4int fNumberOfSteps;
  GXTrack* fpTrack;
  GXeProcessManager* fpProcessManager;
  GXeSteppingManager* fpSteppingManager;

};

#endif
