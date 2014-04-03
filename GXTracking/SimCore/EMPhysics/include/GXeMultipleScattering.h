#ifndef GXeMultipleScattering_h
#define GXeMultipleScattering_h 1

#include "GPGPILSelection.h"
#include "GPForceCondition.h"
#include "GPVParticleChange.h"

#include "GXTrack.h"
#include "GXPhysicsTable.h"

#include "GPRandom.h"
#include "GPMaterial.h"

class GXeMultipleScattering 
{
public:

  FQUALIFIER
  GXeMultipleScattering(curandState* devStates, 
			G4int threadId,
			GXPhysicsTable* lambdaTable);

  FQUALIFIER 
  ~GXeMultipleScattering() {};

  FQUALIFIER
  G4double AlongStepGetPhysicalInteractionLength(GPGPILSelection* selection);

  FQUALIFIER
  G4double PostStepGetPhysicalInteractionLength(GXTrack* track,
                                                G4double previousStepSize,
						GPForceCondition* condition);

  FQUALIFIER
  GPVParticleChange& AlongStepDoIt(GXTrack* track, 
				   GPMaterial* material); 
  FQUALIFIER
  GPVParticleChange& PostStepDoIt(GXTrack* track);

private:

  //curand
  int fThreadId;
  curandState* fDevStates;

  //physics tables
  GXPhysicsTable* theLambdaTable;

  //G4VMultipleScattering
  GPVParticleChange           fParticleChange;

};

#endif
