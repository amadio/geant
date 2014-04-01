#ifndef GXeIonisation_H
#define GXeIonisation_H 1

#include "GPTypeDef.h"
#include "GPGPILSelection.h"
#include "GPForceCondition.h"
#include "GPVParticleChange.h"

#include "GXTrack.h"
#include "GXPhysicsTable.h"

#include "GPRandom.h"
#include "GPMaterial.h"

class GXeIonisation
{
public:

  FQUALIFIER 
  GXeIonisation(curandState* devStates, 
		G4int threadId, 
		GXPhysicsTable* lamdaTable,
		GXPhysicsTable* rangeTable, 
		GXPhysicsTable* dedxTable, 
		GXPhysicsTable* invrTable);

  FQUALIFIER 
  ~GXeIonisation() {};

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

  //random states
  int fThreadId;
  curandState* fDevStates;

  //physics tables
  GXPhysicsTable* theLambdaTable;
  GXPhysicsTable* theRangeTable;
  GXPhysicsTable* theDEDXTable;
  GXPhysicsTable* theInverseRangeTable;

  GPVParticleChange fParticleChange;

};

#endif

