#ifndef GXeBremsstrahlung_H
#define GXeBremsstrahlung_H 1

#include "GPTypeDef.h"
#include "GPGPILSelection.h"
#include "GPForceCondition.h"
#include "GPVParticleChange.h"

#include "GXTrack.h"
#include "GXPhysicsTable.h"

#include "GPMaterial.h"
#include "GPRandom.h"

class GXeBremsstrahlung
{
public:

  FQUALIFIER 
  GXeBremsstrahlung(curandState* devStates, 
		    G4int threadId, 
		    GXPhysicsTable* lamdaTable);

  FQUALIFIER 
  ~GXeBremsstrahlung() {};

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
  G4int fThreadId;
  curandState* fDevStates;

  //lambda table
  GXPhysicsTable* theLambdaTable;

  GPVParticleChange fParticleChange;
};

#endif

