#ifndef GPeIonisation_H
#define GPeIonisation_H 1

#include "GPTypeDef.h"
#include "GPConstants.h"
#include "GPThreeVector.h"
#include "GPMaterial.h"
#include "GPElement.h"
#include "GPPhysicsTable.h"
#include "GPGPILSelection.h"
#include "GPForceCondition.h"
#include "GPVParticleChange.h"

#include "GXTrack.h"
#include "GPMollerBhabhaModel.h"

class GPeIonisation
{
public:

  FQUALIFIER GPeIonisation(curandState* devStates, int threadId, 
			   GPPhysicsTable* lamdaTable,
			   GPPhysicsTable* rangeTable, 
			   GPPhysicsTable* dedxTable, 
			   GPPhysicsTable* invrTable);

  FQUALIFIER ~GPeIonisation();

  //G4eIonisation InitialiseEnergyLossProcess
  FQUALIFIER void InitialiseProcess(GPMollerBhabhaModel* model);

  //G4VProcess
  FQUALIFIER void ResetNumberOfInteractionLengthLeft();
  FQUALIFIER void EndTracking();

  //G4VEnergyLossProcess
  FQUALIFIER void StartTracking();
  FQUALIFIER void SetIonisation(G4bool val);
  FQUALIFIER G4bool IsIonisationProcess();

  FQUALIFIER void DefineMaterial();
  FQUALIFIER G4double GetLambdaForScaledEnergy(G4double e);
  FQUALIFIER G4double GetScaledRangeForScaledEnergy(G4double e);
  FQUALIFIER G4double GetDEDXForScaledEnergy(G4double e);
  FQUALIFIER G4double ScaledKinEnergyForLoss(G4double r);

  FQUALIFIER G4double AlongStepGetPhysicalInteractionLength(
				  GPGPILSelection* selection);

  FQUALIFIER G4double PostStepGetPhysicalInteractionLength(
			          GXTrack* track,
				  G4double previousStepSize,
				  GPForceCondition* condition);

  FQUALIFIER GPVParticleChange& AlongStepDoIt(GXTrack* track,
					      GPMaterial* material,
					      G4double stepLength);

  FQUALIFIER GPVParticleChange& PostStepDoIt(GXTrack* track,
					     GPMaterial* material);

private:

  //curand
  int fThreadId;
  curandState* fDevStates;

  //G4eIonisation
  G4bool   isInitialised;
  GPMollerBhabhaModel* currentModel; //Ionisation Model

  //G4VProcess
  G4double theNumberOfInteractionLengthLeft;
  G4double currentInteractionLength;
  G4double theInitialNumberOfInteractionLength;
  GPVParticleChange fParticleChange;

  //G4VEnergyLossProcess
  G4bool isIonisation;
  GPGPILSelection  aGPILSelection;

  GPPhysicsTable* theLambdaTable;

  G4double preStepLambda;
  G4double preStepKinEnergy;
  G4double preStepScaledEnergy;
  G4double mfpKinEnergy;

  G4double fRange;
  G4double massRatio;
  G4double fFactor;
  G4double reduceFactor;
  G4double lowestKinEnergy;

  G4double dRoverRange;
  G4double finalRange;
  G4double minKinEnergy;

  GPPhysicsTable* theRangeTable;
  GPPhysicsTable* theDEDXTable;
  GPPhysicsTable* theInverseRangeTable;


};

#endif

