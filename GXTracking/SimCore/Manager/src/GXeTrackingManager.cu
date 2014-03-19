#include "GXeTrackingManager.h"

FQUALIFIER
GXeTrackingManager::GXeTrackingManager(GXeProcessManager* aProcessManager,
				       GXeSteppingManager* aSteppingManager)
  :  fNumberOfSteps(1),
     fpTrack(0)
{
  fpProcessManager = aProcessManager;
  fpSteppingManager = aSteppingManager;
}

FQUALIFIER
void GXeTrackingManager::SetNumberOfSteps(G4int aValue)
{
  fNumberOfSteps = aValue;
}

FQUALIFIER
void GXeTrackingManager::ProcessGPIL(GXTrack* apValueGXTrack)
{
  // Receiving a GXTrack from a track scheduler or broker, 
  fpTrack = apValueGXTrack;

  // Give SteppingManger the pointer to the track which will be tracked 
  fpSteppingManager->SetGPILStep(fpTrack);

  // Give SteppingManger the maxmimum number of processes 
  fpSteppingManager->GetProcessNumber(fpProcessManager);

  // Process GPIL
  fpSteppingManager->SteppingGPIL();

}

FQUALIFIER
void GXeTrackingManager::ProcessDoIt(GXTrack* apValueGXTrack)
{
  // Receiving a GXTrack from a track scheduler or broker, 
  fpTrack = apValueGXTrack;

  // Give SteppingManger the pointer to the track which will be tracked 
  //  fpSteppingManager->SetDoItStep(fpTrack,apLiasonGXTrack);
  fpSteppingManager->SetDoItStep(fpTrack);

  // Give SteppingManger the maxmimum number of processes 
  fpSteppingManager->GetProcessNumber(fpProcessManager);

  // Track the particle Step-by-Step while it is alive
  fpSteppingManager->SteppingDoIt();
}


FQUALIFIER 
GXTrack* GXeTrackingManager::GetTrack()
{ 
  return fpTrack;
}

FQUALIFIER 
GXeSteppingManager* GXeTrackingManager::GetSteppingManager()
{ 
  return fpSteppingManager; 
}
