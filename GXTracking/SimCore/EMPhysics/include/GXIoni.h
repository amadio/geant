#ifndef GXIoni_h
#define GXIoni_h

#include "GPPhysicsTable.h"
#include "GPForceCondition.h"

#include "GXTrack.h"
#include "GPThreeVector.h"

#include "GPRandom.h"

struct GXIoni {

public:

  FQUALIFIER GXIoni();
  FQUALIFIER ~GXIoni() {}

  FQUALIFIER void Print();
  FQUALIFIER void SetCurandState(curandState* v);
  FQUALIFIER void UseIntegral(bool v);
  FQUALIFIER void InitialiseStep(G4double kineticEnergy);
  FQUALIFIER void SetLambdaTable(GPPhysicsTable* val);
  FQUALIFIER void UseLambdaTable(bool v);
  FQUALIFIER void ResetNumberOfInteractionLengthLeft();
  FQUALIFIER void SubtractNumberOfInteractionLengthLeft(G4double previousStepSize);
  FQUALIFIER G4double ComputeCrossSectionPerElectron(double kineticEnergy, double cutEnergy, double maxEnergy);
  FQUALIFIER G4double ComputeCrossSectionPerAtom(double kineticEnergy, double Z, double cut);
  FQUALIFIER G4double CrossSectionPerVolume(double kineticEnergy, double cutEnergy, double maxEnergy);
  FQUALIFIER G4double GetCurrentLambda(double e);
  FQUALIFIER void ComputeIntegralLambda(G4double e);
  FQUALIFIER G4double PostStepGetPhysicalInteractionLength(double kineticEnergy, double previousStepSize, GPForceCondition* condition);
  FQUALIFIER G4double SampleSecondaries(GXTrack* dp, double tmin, double maxEnergy);
  FQUALIFIER G4double PostStepDoIt(GXTrack* dp);


  curandState* localState;

  bool useLambdaTable;
  bool integral;

  G4double lambdaFactor;
  GPPhysicsTable* theLambdaTable;

  G4double mfpKinEnergy;
  G4double preStepKinEnergy;
  G4double preStepLambda;
  G4double fFactor;

  // G4VProcess
  G4double theNumberOfInteractionLengthLeft;
  G4double currentInteractionLength;

  // G4VEmModel
  bool isElectron;
  bool LPMFlag;
  G4double xsec[nElements];
};

#endif
