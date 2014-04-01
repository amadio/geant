#ifndef GPeMultipleScattering_h
#define GPeMultipleScattering_h 1

//#include "G4VContinuousDiscreteProcess.hh"
//#include "globals.hh"
//#include "G4ParticleChangeForMSC.hh"
//#include "G4Step.hh"
//#include "G4EmModelManager.hh"
//#include "G4VMscModel.hh"
//#include "G4MscStepLimitType.hh"
//class G4ParticleDefinition;
//class G4VEnergyLossProcess;
//class G4LossTableManager;
//class G4SafetyHelper;

#include "GPMaterial.h"
#include "GPPhysicsTable.h"
#include "GPGPILSelection.h"
#include "GPForceCondition.h"
#include "GPVParticleChange.h"

#include "GXTrack.h"
#include "GPUrbanMscModel95.h"

class GPeMultipleScattering 
{
public:

  FQUALIFIER
  GPeMultipleScattering(curandState* devStates, 
			int threadId,
			GPPhysicsTable* lambdaTable);

  FQUALIFIER
  ~GPeMultipleScattering();

  //G4eMultipleScattering
  FQUALIFIER void InitialiseProcess(GPUrbanMscModel95* model);

  // The function overloads the corresponding function of the base
  // class.It limits the step near to boundaries only
  // and invokes the method GetMscContinuousStepLimit at every step.

  FQUALIFIER void StartTracking();
  FQUALIFIER void EndTracking();

  FQUALIFIER
  G4double AlongStepGetPhysicalInteractionLength(GPMaterial* material,
						 G4double kineticEnergy,
						 G4double  currentMinimalStep,
						 GPGPILSelection* selection);

  // The function overloads the corresponding function of the base
  // class.
  FQUALIFIER
  G4double PostStepGetPhysicalInteractionLength(
                                            GXTrack* track,
					    GPForceCondition* condition);

  // Along step actions
  FQUALIFIER
  GPVParticleChange& AlongStepDoIt(GPMaterial* material, GXTrack* track);

  // Post step actions
  FQUALIFIER
  GPVParticleChange& PostStepDoIt(GXTrack* track);

  FQUALIFIER
  G4double GetTransportMeanFreePath(G4double ekin);

private:

  //curand
  int fThreadId;
  curandState* fDevStates;

  //G4VProcess
  G4double theNumberOfInteractionLengthLeft;
  G4double currentInteractionLength;
  G4double theInitialNumberOfInteractionLength;

  //G4VMultipleScattering

  G4double                    geomMin;
  G4double                    lowestKinEnergy;
  GPVParticleChange           fParticleChange;

  G4double                    physStepLimit;
  G4double                    tPathLength;
  G4double                    gPathLength;

  GPThreeVector               fNewPosition;
  G4bool                      fPositionChanged;
  G4bool                      isActive;

  //G4eMultipleScattering
  G4bool   isInitialized;
  GPUrbanMscModel95* currentModel; 
  GPPhysicsTable* theLambdaTable;

};


#endif
