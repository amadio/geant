#ifndef GXeBremsstrahlung_H
#define GXeBremsstrahlung_H 1

#include "GPTypeDef.h"
#include "GPRandom.h"
#include "GPGPILSelection.h"
#include "GPForceCondition.h"
#include "GPVParticleChange.h"

#include "GXTrack.h"
#include "GXPhysicsTable.h"

#include "GXSeltzerBergerModel.h"

class GXeBremsstrahlung
{
public:

  FQUALIFIER 
  GXeBremsstrahlung(curandState* devStates, 
		    G4int threadId, 
		    GPPhysicsTable* lamdaTable);

  FQUALIFIER 
  GXeBremsstrahlung(curandState* devStates, 
		    G4int threadId, 
		    G4double* xsecTable);
  
  FQUALIFIER 
  ~GXeBremsstrahlung() {};

  FQUALIFIER 
  void InitialiseProcess(GXSeltzerBergerModel* model);

  FQUALIFIER 
  void ResetNumberOfInteractionLengthLeft_G4();

  FQUALIFIER 
  void StartTracking_G4();

  FQUALIFIER 
  void StartTracking_s(GXTrack* track);

  FQUALIFIER 
  void StartTracking_v(G4int ntracks,
		       GXTrack* track);

  FQUALIFIER 
  void StartTracking_soa(G4int ntracks,
			 GXSoATrack track);

  FQUALIFIER 
  void DefineMaterial(GPMaterial* material);
  
  FQUALIFIER 
  G4double AlongStepGPIL(GPGPILSelection* selection);
  
  FQUALIFIER 
  G4double PostStepGPIL_G4(GXTrack* track,
			   GPMaterial* material,
			   G4double previousStepSize,
			   GPForceCondition* condition);

  FQUALIFIER 
  G4double PostStepGPIL_s(GXTrack* track,
			  GPMaterial* material,
			  G4double previousStepSize,
			  GPForceCondition* condition);
  
  FQUALIFIER 
  G4double PostStepGPIL_t(GXTrack* track,
			  GPMaterial* material,
			  G4double previousStepSize,
			  GPForceCondition* condition);

  FQUALIFIER 
  void PostStepGPIL_v(G4int ntracks,
		      GXTrack* track,
		      GPMaterial* material,
		      GPForceCondition* condition);

  FQUALIFIER 
  void PostStepGPIL_soa(G4int ntracks,
			GXSoATrack track,
			GPMaterial* material,
			GPForceCondition* condition);

  FQUALIFIER 
  GPVParticleChange& AlongStepDoIt(GXTrack* track, 
				   GPMaterial* material);

  FQUALIFIER 
  GPVParticleChange& PostStepDoIt_G4(GXTrack* track,
				     GPMaterial* material);

  FQUALIFIER 
  GPVParticleChange& PostStepDoIt_s(GXTrack* track,
				    GPMaterial* material);

  FQUALIFIER 
  GPVParticleChange& PostStepDoIt_t(GXTrack* track,
				    GPMaterial* material);

  FQUALIFIER 
  void PostStepDoIt_v(G4int ntracks,
		      GXTrack* track,
		      GXTrack* secTracks,
		      G4int* stackSize,
		      GPMaterial* material);

  FQUALIFIER 
  G4double GetLambdaForScaledEnergy_G4(G4double e);

  FQUALIFIER 
  G4double GetLambda(G4double energy);

  FQUALIFIER 
  G4double GetLambda_s(GXTrack* track);

  FQUALIFIER 
  G4double GetLambda_t(GXTrack* track);

  FQUALIFIER 
  void GetLambda_v(G4int ntracks,
		   GXTrack* track);

  FQUALIFIER 
  void GetLambda_soa(G4int ntracks,
		     GXSoATrack track);

private:

  //curand
  G4int fThreadId;
  curandState* fDevStates;

  //tabulated cross section data
  G4double* fXsecTable;

  //Geant4 lambda table and variables from G4VProcess
  GPPhysicsTable* fLambdaTable;
  G4double fPreStepLambda;

  //G4eBremsstrahlung Model
  G4bool isInitialised;
  GPVParticleChange fParticleChange;
  GXSeltzerBergerModel* currentModel; 

  //G4VProcess
  G4double theNumberOfInteractionLengthLeft;
  G4double currentInteractionLength;
  G4double theInitialNumberOfInteractionLength;
  G4double preStepLambda;
  G4double preStepKinEnergy;
  G4double preStepScaledEnergy;
  G4double mfpKinEnergy;
  G4double massRatio;
  G4double fFactor;
  G4double lowestKinEnergy;

};

#endif

