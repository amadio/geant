#ifndef GPProcessManager_H
#define GPProcessManager_H 1

#include "GPTypeDef.h"
#include "GPConstants.h"
#include "GPPhotonProcess.h"
#include "GPeBremsstrahlung.h"
#include "GPeIonisation.h"
#include "GPeMultipleScattering.h"

class GPProcessManager 
{
public:
 
  FQUALIFIER GPProcessManager();
  FQUALIFIER ~GPProcessManager();

  FQUALIFIER void AddPhotonProcess(GPPhotonProcess& process);
  FQUALIFIER G4int GetNumberOfPhotonProcess();
  FQUALIFIER GPPhotonProcess* GetPhotonProcessVector();

  FQUALIFIER void AddElectronProcess(GPeBremsstrahlung* aBrem,
				     GPeIonisation* aIoin,
				     GPeMultipleScattering* aMsc);
  FQUALIFIER G4int GetNumberOfElectronProcess();
  FQUALIFIER GPeBremsstrahlung* GetBremsstrahlungProcess();
  FQUALIFIER GPeIonisation* GetIonisationProcess();
  FQUALIFIER GPeMultipleScattering* GetMultipleScatteringProcess();

  FQUALIFIER void StartTracking(GXTrack* aTrack=0);
  FQUALIFIER void EndTracking();

private:
  G4int theNumElectronProcesses;
  G4int theNumPhotonProcesses;

  GPPhotonProcess thePhotonProcVector[nPhotonProcessMax];
  GPeBremsstrahlung* theBremsstrahlungProcess;
  GPeIonisation* theIonisationProcess;
  GPeMultipleScattering* theMultipleScatteringProcess;
};

#endif
