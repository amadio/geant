#include "GXeBremsstrahlung.h"
#include "GPRandom.h"

FQUALIFIER
GXeBremsstrahlung::GXeBremsstrahlung(curandState* devStates,
				     G4int threadId,
				     GXPhysicsTable* lamdaTable)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  //lambda table
  theLambdaTable = lamdaTable;

}

FQUALIFIER
G4double GXeBremsstrahlung::AlongStepGetPhysicalInteractionLength(
					        GPGPILSelection* selection)
{
  // condition is set to "Not Forced" 
  *selection = NotCandidateForSelection;

  // get GPIL
  G4double x = DBL_MAX;
  return x;
}

FQUALIFIER
G4double GXeBremsstrahlung::PostStepGetPhysicalInteractionLength(GXTrack* track,
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
GPVParticleChange& GXeBremsstrahlung::AlongStepDoIt(GXTrack* track, 
						    GPMaterial* material) 
{
  fParticleChange.ParticleChangeForLoss_InitializeForAlongStep(track);
  return fParticleChange;
}

FQUALIFIER
GPVParticleChange& GXeBremsstrahlung::PostStepDoIt(GXTrack* track)
{
  fParticleChange.ParticleChangeForLoss_InitializeForPostStep(track);
  return fParticleChange;
}
