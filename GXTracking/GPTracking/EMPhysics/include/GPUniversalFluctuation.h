#ifndef GPUniversalFluctuation_h
#define GPUniversalFluctuation_h 1

#include "GPMaterial.h"

#include "GPRandom.h"

class GPUniversalFluctuation 
{

public:

  FQUALIFIER GPUniversalFluctuation(curandState* devStates, 
				    int threadId);

  FQUALIFIER ~GPUniversalFluctuation();

  FQUALIFIER G4double SampleFluctuations(GPMaterial* material,
					 G4double kineticEnergy,
					 G4double& tmax,
					 G4double& length,
					 G4double&averageLoss);

  FQUALIFIER G4double Dispersion(GPMaterial* material,
				 G4double kineticEnergy,
				 G4double& tmax,
				 G4double& length);

  FQUALIFIER void InitialiseMe();

  // Initialisation prestep
  FQUALIFIER void SetParticleAndCharge(G4double q2);

private:

  int fThreadId;
  curandState* fDevStates;

  // hide assignment operator
  //  GPUniversalFluctuation & operator=(const  GPUniversalFluctuation &right);
  //  GPUniversalFluctuation(const  GPUniversalFluctuation&);

  //  const G4ParticleDefinition* particle;
  GPMaterial* lastMaterial;

  G4double particleMass;
  G4double chargeSquare;

  // data members to speed up the fluctuation calculation
  G4double ipotFluct;
  G4double electronDensity;
  
  G4double f1Fluct;
  G4double f2Fluct;
  G4double e1Fluct;
  G4double e2Fluct;
  G4double e1LogFluct;
  G4double e2LogFluct;
  G4double ipotLogFluct;
  G4double e0;
  G4double esmall;

  G4double e1,e2;

  G4double minNumberInteractionsBohr;
  G4double theBohrBeta2;
  G4double minLoss;
  G4double nmaxCont;
  G4double rate,fw;

};

#endif

