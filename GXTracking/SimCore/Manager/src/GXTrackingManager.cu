#include "GXTrackingManager.h"

FQUALIFIER
GXTrackingManager::GXTrackingManager(GXProcessManager* aProcessManager,
				     GXSteppingManager* aSteppingManager)
  :  fNumberOfSteps(1),
     EventIsAborted(false)
{
  fpProcessManager = aProcessManager;
  fpSteppingManager = aSteppingManager;
}

FQUALIFIER
void GXTrackingManager::SetNumberOfSteps(G4int aValue)
{
  fNumberOfSteps = aValue;
}

FQUALIFIER
void GXTrackingManager::ProcessOneTrack(GXTrack* apValueG4Track)
{
  // Receiving a GXTrack from a track scheduler or broker, 
  // this funciton has the responsibility to trace the track till it stops.
  fpTrack = apValueG4Track;
  EventIsAborted = false;

  // @@@G4FWP handle Secondaries via a container on the global memory
  // Clear 2ndary particle vector
  //  size_t itr;
  //  for(itr=0;itr<GimmeSecondaries()->size();itr++){ 
  //     delete (*GimmeSecondaries())[itr];
  //  }
  //  GimmeSecondaries()->clear();  
   
  // Give SteppingManger the pointer to the track which will be tracked 
  fpSteppingManager->SetInitialStep(fpTrack);

  // Pre tracking user intervention process.

  // Give SteppingManger the maxmimum number of processes 
  fpSteppingManager->GetProcessNumber(fpProcessManager);

  // Give track the pointer to the Step 
  // @@@G4FWP GXTrack doesn't own a GPStep - set Stet via GXSteppingManager
  //  fpTrack->SetStep(fpSteppingManager->GetStep());

  // Inform beginning of tracking to physics processes 
  fpProcessManager->StartTracking(fpTrack);

  // Track the particle Step-by-Step while it is alive

  G4int  stepCounter = 0;
  while( ((fpTrack->status == fAlive) ||
	  (fpTrack->status == fStopButAlive)) && 
         stepCounter < fNumberOfSteps ){

    //    fpTrack->IncrementCurrentStepNumber();
    ++stepCounter;
    fpSteppingManager->Stepping();

    if(EventIsAborted) {
      fpTrack->status = fKillTrackAndSecondaries;
    }
  }

  // Inform end of tracking to physics processes 
  fpProcessManager->EndTracking();

}

FQUALIFIER
void GXTrackingManager::ProcessGPIL(GXTrack* apValueGXTrack)
{
  // Receiving a GXTrack from a track scheduler or broker, 
  fpTrack = apValueGXTrack;

  // Give SteppingManger the pointer to the track which will be tracked 
  //  fpSteppingManager->SetGPILStep(fpTrack,apLiasonGXTrack);
  fpSteppingManager->SetGPILStep(fpTrack);

  // Give SteppingManger the maxmimum number of processes 
  fpSteppingManager->GetProcessNumber(fpProcessManager);

  // Inform beginning of tracking to physics processes 
  fpProcessManager->StartTracking(fpTrack);

  // Process GPIL
  fpSteppingManager->SteppingGPIL();

}

FQUALIFIER
void GXTrackingManager::ProcessDoIt(GXTrack* apValueGXTrack)
{
  // Receiving a GXTrack from a track scheduler or broker, 
  fpTrack = apValueGXTrack;
  EventIsAborted = false;

  // Give SteppingManger the pointer to the track which will be tracked 
  //  fpSteppingManager->SetDoItStep(fpTrack,apLiasonGXTrack);
  fpSteppingManager->SetDoItStep(fpTrack);

  // Give SteppingManger the maxmimum number of processes 
  fpSteppingManager->GetProcessNumber(fpProcessManager);

  // Track the particle Step-by-Step while it is alive

  fpSteppingManager->SteppingDoIt();

  if(EventIsAborted) {
    fpTrack->status = fKillTrackAndSecondaries;
  }

  // Inform end of tracking to physics processes 
  fpProcessManager->EndTracking();
}

FQUALIFIER
void GXTrackingManager::EventAborted()
{
  fpTrack->status = fKillTrackAndSecondaries;
  EventIsAborted = true;
}

FQUALIFIER 
GXTrack* GXTrackingManager::GetTrack()
{ 
  return fpTrack;
}

FQUALIFIER 
GXSteppingManager* GXTrackingManager::GetSteppingManager()
{ 
  return fpSteppingManager; 
}
