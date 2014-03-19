#include "GXeProcessManager.h"
#include <assert.h>

FQUALIFIER 
GXeProcessManager::GXeProcessManager()
{
  theNumElectronProcesses = 0;

  theBremsstrahlungProcess = 0;   
  theIonisationProcess = 0;
  theMultipleScatteringProcess = 0;
}

FQUALIFIER 
GXeProcessManager::~GXeProcessManager()
{
  ;
}

FQUALIFIER 
void GXeProcessManager::AddElectronProcess(GXeBremsstrahlung* aBrem,
					   GXeIonisation* aIoni,
					   GXeMultipleScattering* aMsc)
{ 
  if(aBrem) {
    theBremsstrahlungProcess = aBrem;
    ++theNumElectronProcesses;
  }
  if(aIoni) {
    theIonisationProcess = aIoni;
    ++theNumElectronProcesses;
  }
  if(aMsc) {
    theMultipleScatteringProcess = aMsc;
    ++theNumElectronProcesses;
  }
}

FQUALIFIER 
GXeBremsstrahlung* GXeProcessManager::GetBremsstrahlungProcess()
{
  return theBremsstrahlungProcess;
}

FQUALIFIER 
GXeIonisation* GXeProcessManager::GetIonisationProcess()
{
  return theIonisationProcess;
}

FQUALIFIER 
GXeMultipleScattering* GXeProcessManager::GetMultipleScatteringProcess()
{
  return theMultipleScatteringProcess;
}

FQUALIFIER 
G4int GXeProcessManager::GetNumberOfElectronProcess()
{
  return theNumElectronProcesses;
}

