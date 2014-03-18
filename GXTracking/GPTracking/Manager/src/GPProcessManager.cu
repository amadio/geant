#include "GPProcessManager.h"
#include <assert.h>

FQUALIFIER 
GPProcessManager::GPProcessManager()
{
  theNumElectronProcesses = 0;
  theNumPhotonProcesses = 0;

  theBremsstrahlungProcess = 0;   
  theIonisationProcess = 0;
  theMultipleScatteringProcess = 0;
}

FQUALIFIER 
GPProcessManager::~GPProcessManager()
{
}

FQUALIFIER 
void GPProcessManager::AddPhotonProcess(GPPhotonProcess& aProcess)
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
void GPProcessManager::AddElectronProcess(GPeBremsstrahlung* aBrem,
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
GPeBremsstrahlung* GPProcessManager::GetBremsstrahlungProcess()
{
  return theBremsstrahlungProcess;
}

FQUALIFIER 
GPeIonisation* GPProcessManager::GetIonisationProcess()
{
  return theIonisationProcess;
}

FQUALIFIER 
GPeMultipleScattering* GPProcessManager::GetMultipleScatteringProcess()
{
  return theMultipleScatteringProcess;
}

FQUALIFIER 
G4int GPProcessManager::GetNumberOfElectronProcess()
{
  return theNumElectronProcesses;
}

FQUALIFIER 
G4int GPProcessManager::GetNumberOfPhotonProcess()
{
  return theNumPhotonProcesses;
}

FQUALIFIER
GPPhotonProcess* GPProcessManager::GetPhotonProcessVector() 
{
  return &(thePhotonProcVector[0]);
}

FQUALIFIER 
void GPProcessManager::StartTracking(GXTrack* aTrack)
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
void GPProcessManager::EndTracking()
{
  for (G4int idx = 0; idx < theNumPhotonProcesses ; idx++){
    thePhotonProcVector[idx].EndTracking();
  }
  if(theBremsstrahlungProcess) theBremsstrahlungProcess->EndTracking();   
  if(theIonisationProcess) theIonisationProcess->EndTracking();
  if(theMultipleScatteringProcess) theMultipleScatteringProcess->EndTracking();
}
