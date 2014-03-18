#include "GXeIonisation.h"
#include "GPRandom.h"

FQUALIFIER
GXeIonisation::GXeIonisation(curandState* devStates,
			     G4int threadId,
			     GXPhysicsTable* lamdaTable, 
			     GXPhysicsTable* rangeTable, 
			     GXPhysicsTable* dedxTable, 
			     GXPhysicsTable* invrTable)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  //physics tables
  theLambdaTable = lamdaTable;
  theRangeTable = rangeTable;
  theDEDXTable = dedxTable;
  theInverseRangeTable = invrTable;
}

FQUALIFIER
G4double GXeIonisation::AlongStepGetPhysicalInteractionLength(
					    GPGPILSelection* selection)
{
  // condition is set to "CandidateForSelection" 
  *selection = CandidateForSelection;

  // get GPIL
  G4double x = DBL_MAX;

  //@@@ add implementation for ionization

  return x;
}

FQUALIFIER
G4double GXeIonisation::PostStepGetPhysicalInteractionLength(GXTrack* track,
					          G4double previousStepSize,
					          GPForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  // get GPIL
  G4double preStepLambda = theLambdaTable->physicsVectors[1].Value(track->E);
  G4double x = -log(GPUniformRand(fDevStates,fThreadId))/preStepLambda;

  return x;
}

FQUALIFIER
GPVParticleChange& GXeIonisation::AlongStepDoIt(GXTrack* track, 
						GPMaterial* material) 
{
  fParticleChange.ParticleChangeForLoss_InitializeForAlongStep(track);
  return fParticleChange;
}

FQUALIFIER
GPVParticleChange& GXeIonisation::PostStepDoIt(GXTrack* track)
{
  fParticleChange.ParticleChangeForLoss_InitializeForPostStep(track);
  return fParticleChange;
}
