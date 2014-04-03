#ifndef GPMollerBhabhaModel_h
#define GPMollerBhabhaModel_h 1

#include "GPMaterial.h"
#include "GPThreeVector.h"
#include "GXTrack.h"

#include "GPVParticleChange.h"

#include "GPRandom.h"

class GPMollerBhabhaModel
{

public:

  FQUALIFIER GPMollerBhabhaModel(curandState* devStates, int threadId);

  FQUALIFIER ~GPMollerBhabhaModel();

  FQUALIFIER G4double MaxSecondaryEnergy(G4double kinEnergy);

  FQUALIFIER void SetElectron(G4bool eflag);

  FQUALIFIER G4double ComputeCrossSectionPerElectron(G4double kineticEnergy,
						     G4double cutEnergy,
						     G4double maxEnergy);
				 
  FQUALIFIER G4double ComputeCrossSectionPerAtom(G4double kineticEnergy,
						 G4double Z,
						 G4double cutEnergy,
						 G4double maxEnergy);
				 				 
  FQUALIFIER G4double CrossSectionPerVolume(GPMaterial* material,
					    G4double kineticEnergy,
					    G4double cutEnergy,
					    G4double maxEnergy);
				 
  FQUALIFIER G4double ComputeDEDXPerVolume(GPMaterial* material,
					   G4double kineticEnergy,
					   G4double cutEnergy);

  FQUALIFIER void SampleSecondaries(GXTrack* track,
				    GPMaterial* material,
				    G4double tmin,
				    G4double maxEnergy);

  FQUALIFIER GXTrack& GetSecondary();

  FQUALIFIER void UpdateTrackMomentum(GXTrack* track,
				      GPThreeVector direction,
				      G4double energy);

  FQUALIFIER void FillSecondary(GXTrack* track, GPThreeVector direction,
				G4double energy, G4double charge);

  //  const G4ParticleDefinition* particle;
  //  G4ParticleDefinition*       theElectron;
  //  G4ParticleChangeForLoss*    fParticleChange;

  //G4VEmModel
  FQUALIFIER void SetParticleChange(GPVParticleChange* p);
  FQUALIFIER G4double HighEnergyLimit();
  FQUALIFIER G4double LowEnergyLimit();
  FQUALIFIER G4double HighEnergyActivationLimit();
  FQUALIFIER G4double LowEnergyActivationLimit();
  FQUALIFIER void SetHighEnergyLimit(G4double val);
  FQUALIFIER void SetLowEnergyLimit(G4double val);
  FQUALIFIER void SetActivationHighEnergyLimit(G4double val);
  FQUALIFIER void SetActivationLowEnergyLimit(G4double val);
  FQUALIFIER G4bool IsActive(G4double kinEnergy);

private:
  //curand                  
  int fThreadId;
  curandState* fDevStates;

  G4bool   isElectron;
  G4double twoln10;
  G4double lowLimit;

  GXTrack theDeltaRay;

  //G4VEmModel
  GPVParticleChange* pParticleChange;
  G4double        highLimit;
  G4double        eMinActive;
  G4double        eMaxActive;

};

#endif
