#ifndef GPSteppingManager_H
#define GPSteppingManager_H 1

#include "GPTypeDef.h"
#include "GPStep.h"
#include "GXTrack.h"
#include "GPMaterial.h"
#include "GPStepStatus.h"
#include "GPForceCondition.h"
#include "GPGPILSelection.h"

#include "GPProcessManager.h"
#include "GPPhotonProcess.h"
#include "GPGeomManager.h"
#include "GPFieldMap.h"
#include "GPNavigator.h"

#include "GXTransportation.h"

class GPSteppingManager 
{
public:

  FQUALIFIER GPSteppingManager(GPGeomManager *geomManager,
			       GPFieldMap *magMap);
  FQUALIFIER ~GPSteppingManager();
  
  FQUALIFIER
    void SetTransportation(GXTransportation* aT,
			   GPNavigator* aNavigator);

  FQUALIFIER GPStepStatus Stepping();
  FQUALIFIER void SetInitialStep(GXTrack* valueTrack);

  FQUALIFIER void SetStep(GPStep* aStep);
  FQUALIFIER GPStep* GetStep();

  FQUALIFIER void SetMaterial(GPMaterial* aMaterial);
  FQUALIFIER GPMaterial* GetMaterial();

  FQUALIFIER void GetProcessNumber(GPProcessManager* aProcessManager);

  FQUALIFIER void SetSecondaryStack(GXTrack *secTracks,
				    G4int *stackSize,
				    G4int *offset);

private:

  FQUALIFIER void DefinePhysicalStepLength();
  FQUALIFIER void InvokeAlongStepDoItProcs();
  FQUALIFIER void InvokePostStepDoItProcs();
  FQUALIFIER void InvokePSDIP(size_t);
  FQUALIFIER void InvokePSDIPeBremsstrahlung();
  FQUALIFIER void InvokePSDIPeIonisation();
  FQUALIFIER void InvokePSDIPeMultipleScattering();
  FQUALIFIER G4double CalculateSafety();

  FQUALIFIER void StoreSecondary(GXTrack* aTrack);
  
private:

  GPStepStatus fStepStatus;
  GXTrack* fTrack;
  GPStep* fStep;
  GPMaterial* fMaterial; //@@@G4FWP replace this by a navigator or touchable
  //  G4TouchableHandle fTouchableHandle;
  GPNavigator* fNavigator;
  GPVPhysicalVolume* fCurrentVolume;

  G4int fN2ndariesPostStepDoIt;

  G4double PhysicalStep;
  G4double GeometricalStep;
  G4double CorrectedStep;
  G4double fPreviousStepSize;
  
  GPPhotonProcess* fPhotonProcessVector;

  GPeBremsstrahlung* fBrem;
  GPeIonisation* fIoni;
  GPeMultipleScattering* fMsc;

  size_t PhotonMAXofPostStepLoops;
  size_t PhotonMAXofAlongStepLoops;
  size_t ElectronMAXofPostStepLoops;
  size_t ElectronMAXofAlongStepLoops;

  //  size_t fAlongStepDoItProcTriggered;
  size_t fPostStepDoItProcTriggered;

  G4double proposedSafety;
  G4double endpointSafety;
  GPThreeVector endpointSafOrigin;

  G4double physIntLength;
  GPForceCondition fCondition;
  GPGPILSelection  fGPILSelection;

  GPVParticleChange fParticleChange;
  GPPhotonProcess* fCurrentProcess;

  //secondary 
  GXTrack *theSecondaryStack;
  G4int *theStackSize;
  G4int *theOffset;

  GXTransportation* fTransportationProcess;
};

#endif
