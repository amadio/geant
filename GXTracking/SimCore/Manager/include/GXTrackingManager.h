#ifndef GXTrackingManager_H
#define GXTrackingManager_H 1

#include "GPTypeDef.h"
#include "GXProcessManager.h"       
#include "GXSteppingManager.h"      
#include "GXTrack.h"                
#include "GPTrackStatus.h"          
#include "GPStepStatus.h"           

class GXTrackingManager 
{
public:

  FQUALIFIER GXTrackingManager(GXProcessManager* aProcessManager,
			       GXSteppingManager* aSteppingManager);

  FQUALIFIER ~GXTrackingManager() {};

  FQUALIFIER GXTrack* GetTrack();
  FQUALIFIER GXSteppingManager* GetSteppingManager();

  FQUALIFIER void SetNumberOfSteps(G4int aValue);
  FQUALIFIER void ProcessOneTrack(GXTrack* apValue);
  FQUALIFIER void ProcessGPIL(GXTrack* apValue);
  FQUALIFIER void ProcessDoIt(GXTrack* apValue);

  FQUALIFIER void EventAborted();

private:

  GXTrack* fpTrack;
  GXProcessManager* fpProcessManager;
  GXSteppingManager* fpSteppingManager;
  G4int fNumberOfSteps;
  G4bool EventIsAborted;

};

#endif
