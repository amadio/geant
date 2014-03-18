#include "GXeMultipleScattering.h"
#include "GPRandom.h"

FQUALIFIER
GXeMultipleScattering::GXeMultipleScattering(curandState* devStates,
					     G4int threadId,
					     GXPhysicsTable* lambdaTable)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  //physics tables
  theLambdaTable = lambdaTable;

}

FQUALIFIER
G4double GXeMultipleScattering::AlongStepGetPhysicalInteractionLength(
						   GPGPILSelection* selection)
{
  // get Step limit proposed by the process
  *selection = NotCandidateForSelection;

  // get GPIL
  G4double x = DBL_MAX;

  //@@@ add implementation for ionization

  return x;
}

FQUALIFIER 
G4double 
GXeMultipleScattering::PostStepGetPhysicalInteractionLength(GXTrack* track, 
                                                G4double previousStepSize,
					        GPForceCondition* condition)
{
  // condition is set to "Forced" 
  *condition = Forced;

  // get GPIL
  return DBL_MAX;
}

FQUALIFIER
GPVParticleChange& GXeMultipleScattering::AlongStepDoIt(GXTrack* track, 
							GPMaterial* material) 
{
  GPThreeVector direction = 
    GPThreeVector_unit(GPThreeVector_create(track->px,track->py,track->pz));
  GPThreeVector fNewPosition = GPThreeVector_create(track->x,track->y,track->z);

  fParticleChange.SetProposedMomentumDirection(direction);
  fParticleChange.SetProposedPosition(fNewPosition);

  //@@@add

  return fParticleChange;

}

FQUALIFIER
GPVParticleChange& GXeMultipleScattering::PostStepDoIt(GXTrack* track)
{  
  fParticleChange.ParticleChangeForMSC_Initialize(track);  
  return fParticleChange;
}
