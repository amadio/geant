#ifndef GPBrem_h
#define GPBrem_h

#include "GPPhysicsTable.h"
#include "GPForceCondition.h"

#include "GPTrack.h"
#include "GPThreeVector.h"

struct GPBrem {

public:

  FQUALIFIER GPBrem();
  FQUALIFIER ~GPBrem() {}

  FQUALIFIER void Print();
  FQUALIFIER void SetCurandStates(curandState* v);
  FQUALIFIER void UseIntegral(bool v);
  FQUALIFIER void InitialiseStep(G4double kineticEnergy);
  FQUALIFIER void SetLambdaTable(GPPhysicsTable* val);
  FQUALIFIER void UseLambdaTable(bool v);
  FQUALIFIER void ResetNumberOfInteractionLengthLeft();
  FQUALIFIER void SubtractNumberOfInteractionLengthLeft(G4double previousStepSize);
  FQUALIFIER G4double SupressionFunction(double kineticEnergy, double gammaEnergy, bool LPMFlag);
  FQUALIFIER G4double PositronCorrFactorSigma(double Z, double kineticEnergy, double cut);
  FQUALIFIER G4double ComputeCrossSectionPerAtom(double kineticEnergy, double Z, double cut);
  FQUALIFIER G4double CrossSectionPerVolume(double kineticEnergy, double cutEnergy, double maxEnergy);
  FQUALIFIER G4double GetCurrentLambda(double e);
  FQUALIFIER void ComputeIntegralLambda(G4double e);
  FQUALIFIER G4double PostStepGetPhysicalInteractionLength(double kineticEnergy, double previousStepSize, GPForceCondition* condition);
  FQUALIFIER G4double ScreenFunction1(double ScreenVariable);
  FQUALIFIER G4double ScreenFunction2(double ScreenVariable);
  FQUALIFIER G4double PolarAngle(double initial_energy, double final_energy, int Z);
  FQUALIFIER G4double GPVEmModel_CrossSectionPerVolume(double ekin, double emin, double emax);
  FQUALIFIER int SelectRandomAtom(G4double kinEnergy, double tcut, double tmax);
  FQUALIFIER G4double SampleSecondaries(GPTrack* dp, double tmin, double maxEnergy);
  FQUALIFIER G4double PostStepDoIt(GPTrack* dp);


  int threadId;
  curandState* devStates;

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
