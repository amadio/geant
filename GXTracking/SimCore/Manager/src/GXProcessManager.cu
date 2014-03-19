#include "GXProcessManager.h"
#include <assert.h>

FQUALIFIER 
GXProcessManager::GXProcessManager()
{
  theNumElectronProcesses = 0;
  theNumPhotonProcesses = 0;

  theBremsstrahlungProcess = 0;   
  theIonisationProcess = 0;
  theMultipleScatteringProcess = 0;
}

FQUALIFIER 
GXProcessManager::~GXProcessManager()
{
  ;
}

FQUALIFIER 
void GXProcessManager::AddPhotonProcess(GPPhotonProcess& aProcess)
{ 
  if(theNumPhotonProcesses < nPhotonProcessMax) {
    thePhotonProcVector[theNumPhotonProcesses] = aProcess;
    theNumPhotonProcesses += 1;
  }
  else {
    assert(0);
  }
}

FQUALIFIER 
void GXProcessManager::AddElectronProcess(GPeBremsstrahlung* aBrem,
					  GPeIonisation* aIoni,
					  GPeMultipleScattering* aMsc)
{ 
  if(aBrem) {
    theBremsstrahlungProcess = aBrem;
    theNumElectronProcesses += 1;
  }
  if(aIoni) {
    theIonisationProcess = aIoni;
    theNumElectronProcesses += 1;
  }
  if(aMsc) {
    theMultipleScatteringProcess = aMsc;
    theNumElectronProcesses += 1;
  }
}

FQUALIFIER 
GPeBremsstrahlung* GXProcessManager::GetBremsstrahlungProcess()
{
  return theBremsstrahlungProcess;
}

FQUALIFIER 
GPeIonisation* GXProcessManager::GetIonisationProcess()
{
  return theIonisationProcess;
}

FQUALIFIER 
GPeMultipleScattering* GXProcessManager::GetMultipleScatteringProcess()
{
  return theMultipleScatteringProcess;
}

FQUALIFIER 
G4int GXProcessManager::GetNumberOfElectronProcess()
{
  return theNumElectronProcesses;
}

FQUALIFIER 
G4int GXProcessManager::GetNumberOfPhotonProcess()
{
  return theNumPhotonProcesses;
}

FQUALIFIER
GPPhotonProcess* GXProcessManager::GetPhotonProcessVector() 
{
  return &(thePhotonProcVector[0]);
}

FQUALIFIER 
void GXProcessManager::StartTracking(GXTrack* aTrack)
{  
  for (G4int idx = 0; idx < theNumPhotonProcesses ; idx++){
    thePhotonProcVector[idx].StartTracking(aTrack);
  }

  if(theBremsstrahlungProcess) theBremsstrahlungProcess->StartTracking();   
  if(theIonisationProcess) theIonisationProcess->StartTracking();
  if(theMultipleScatteringProcess) 
    theMultipleScatteringProcess->StartTracking();
}

FQUALIFIER 
void GXProcessManager::EndTracking()
{
  for (G4int idx = 0; idx < theNumPhotonProcesses ; idx++){
    thePhotonProcVector[idx].EndTracking();
  }
  if(theBremsstrahlungProcess) theBremsstrahlungProcess->EndTracking();   
  if(theIonisationProcess) theIonisationProcess->EndTracking();
  if(theMultipleScatteringProcess) theMultipleScatteringProcess->EndTracking();
}
