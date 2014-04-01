#ifndef GXeProcessManager_H
#define GXeProcessManager_H 1

#include "GPTypeDef.h"

#include "GXeBremsstrahlung.h"
#include "GXeIonisation.h"
#include "GXeMultipleScattering.h"

class GXeProcessManager 
{
public:
 
  FQUALIFIER GXeProcessManager();
  FQUALIFIER ~GXeProcessManager();

  FQUALIFIER void AddElectronProcess(GXeBremsstrahlung* aBrem,
				     GXeIonisation* aIoin,
				     GXeMultipleScattering* aMsc);

  FQUALIFIER G4int GetNumberOfElectronProcess();

  FQUALIFIER GXeBremsstrahlung* GetBremsstrahlungProcess();
  FQUALIFIER GXeIonisation* GetIonisationProcess();
  FQUALIFIER GXeMultipleScattering* GetMultipleScatteringProcess();

private:
  G4int theNumElectronProcesses;

  GXeBremsstrahlung* theBremsstrahlungProcess;
  GXeIonisation* theIonisationProcess;
  GXeMultipleScattering* theMultipleScatteringProcess;
};

#endif
