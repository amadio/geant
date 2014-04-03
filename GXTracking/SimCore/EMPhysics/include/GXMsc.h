#ifndef GXMsc_h
#define GXMsc_h

#include "GPPhysicsTable.h"
#include "GPForceCondition.h"
#include "GPGPILSelection.h"
#include "GPMscStepLimitType.h"
#include "GPStepStatus.h"

#include "GXTrack.h"
#include "GPThreeVector.h"

#include "GPRandom.h"

struct GXMsc {

public:

  FQUALIFIER GXMsc();
  FQUALIFIER ~GXMsc() {}

  FQUALIFIER void Print();
  FQUALIFIER void SetCurandState(curandState* v);
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
  FQUALIFIER G4double AlongStepGetPhysicalInteractionLength(GXTrack* dp, double currentMinimalStep, GPGPILSelection* selection);
  FQUALIFIER G4double ComputeTruePathLengthLimit(GXTrack* dp, double currentMinimalStep);
  FQUALIFIER G4double ComputeGeomPathLength();
  FQUALIFIER G4double GetRange(double e);
  FQUALIFIER G4double PostStepGetPhysicalInteractionLength(double kineticEnergy, double previousStepSize, GPForceCondition* condition);
  FQUALIFIER G4double SampleScattering(GXTrack* dp);
  FQUALIFIER G4double AlongStepDoIt(GXTrack* dp);
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

  // G4VMscModel
  GPMscStepLimitType steppingAlgorithm;
  G4double skin;
  G4double facrange;
  G4double facgeom;
  G4double facsafety;
  G4double dtrl;

  // G4UrbanMscModel95
  bool inside;
  bool insideskin;
  bool samplez;
  G4double currentKinEnergy;
  G4double tPathLength;
  G4double currentRange;
  G4double lambda0;
  G4double presafety;
  G4double smallstep;
  G4double tlimitminfix;
  G4double rangeinit;
  G4double tlimit;
  G4double tlimitmin;
  G4double stepmin;
  G4double skindepth;
  G4double tgeom;
  G4double geombig;
  G4double geommin;
  G4double geomlimit;
  G4double fr;
  G4double masslimite;
  G4double mass;
  G4double lambdalimit;
  G4double par1, par2, par3;
  G4double zPathLength;
  G4double taulim;
  G4double tausmall;
  G4double lambdaeff;
  G4double third;
};

#endif
