#ifndef GPTrackingManager_H
#define GPTrackingManager_H 1

#include "GPTypeDef.h"
#include "GPProcessManager.h"       
#include "GPSteppingManager.h"      
#include "GXTrack.h"                
#include "GPTrackStatus.h"          
#include "GPStepStatus.h"           

class GPTrackingManager 
{
public:

  FQUALIFIER GPTrackingManager(GPProcessManager* aProcessManager,
			       GPSteppingManager* aSteppingManager);

  FQUALIFIER ~GPTrackingManager();

  FQUALIFIER GXTrack* GetTrack();
  FQUALIFIER GPSteppingManager* GetSteppingManager();

  FQUALIFIER void SetNumberOfSteps(G4int aValue);
  FQUALIFIER void ProcessOneTrack(GXTrack* apValue);

  FQUALIFIER void EventAborted();

private:

  GXTrack* fpTrack;
  GPProcessManager* fpProcessManager;
  GPSteppingManager* fpSteppingManager;
  G4int fNumberOfSteps;
  G4bool EventIsAborted;

};

#endif
