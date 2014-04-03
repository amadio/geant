#ifndef GPSeltzerBergerRelModel_H
#define GPSeltzerBergerRelModel_H 1

#include "GPTypeDef.h"
#include "GPConstants.h"
#include "GPThreeVector.h"

#include "GPPhysicsTable.h"
#include "GPPhysics2DVector.h"
#include "GPVParticleChange.h"

#include "GXTrack.h"
#include "GPMaterial.h"
#include "GPElement.h"

#include "GPRandom.h"

class GPSeltzerBergerRelModel
{
public:

  FQUALIFIER GPSeltzerBergerRelModel(curandState* devStates,
				     int threadId,
				     GPPhysics2DVector* data);

  FQUALIFIER ~GPSeltzerBergerRelModel();

  FQUALIFIER GPElement* SelectRandomAtom(GPMaterial* material,
					 G4double kinEnergy,
					 G4double tcut,
					 G4double tmax);


  FQUALIFIER G4double CrossSectionPerVolume(GPMaterial* material,
					    G4double ekin,
					    G4double emin,
					    G4double emax);

  FQUALIFIER void SetupForMaterial(GPMaterial* mat, 
				   G4double kineticEnergy);


  FQUALIFIER G4double ComputeCrossSectionPerAtom(G4double kineticEnergy, 
						 G4double Z,
						 G4double cutEnergy, 
						 G4double maxEnergy);

  FQUALIFIER void SetCurrentElement(G4double Z);

  FQUALIFIER G4double ComputeXSectionPerAtom(G4double cut);

  FQUALIFIER G4double ComputeRelDXSectionPerAtom(G4double gammaEnergy);

  FQUALIFIER G4double ComputeDXSectionPerAtom(G4double gammaEnergy);

  FQUALIFIER G4double ComputeDXSectionPerAtom_SeltzerBergerModel(
					      G4double gammaEnergy);

  FQUALIFIER G4double ComputeDXSectionPerAtom_eBremsstrahlungRelModel(
					      G4double gammaEnergy);

  FQUALIFIER void CalcLPMFunctions(G4double k);

  FQUALIFIER G4double Phi1(G4double gg);

  FQUALIFIER G4double Phi1M2(G4double gg);

  FQUALIFIER G4double Psi1(G4double eps);

  FQUALIFIER G4double Psi1M2(G4double eps);

  FQUALIFIER GPThreeVector DipBustGenerator_SampleDirection(G4double eTkin, 
			  	               GPThreeVector refDirection);

  FQUALIFIER void SampleSecondaries(GXTrack* track,
				    GPMaterial* material,
				    G4double cutEnergy,
				    G4double maxEnergy);

  FQUALIFIER void SampleSecondaries_SeltzerBergerModel(GXTrack* track,
                                                 GPMaterial* material,
                                                 G4double cutEnergy,
						 G4double maxEnergy);

  FQUALIFIER void SampleSecondaries_eBremsstrahlungRelModel(GXTrack* track,
                                                 GPMaterial* material,
                                                 G4double cutEnergy,
						 G4double maxEnergy);

  FQUALIFIER GXTrack& GetSecondary();

  FQUALIFIER void FillSecondary(GXTrack* track, GPThreeVector direction,
                                G4double energy, G4double charge);

  FQUALIFIER void UpdateTrackMomentum(GXTrack* track, GPThreeVector direction,
                                      G4double energy);

  //G4VEmModel
  FQUALIFIER void SetParticleChange(GPVParticleChange* p);
  FQUALIFIER void SetEnergyLimitModes(G4double energyLimit);
  FQUALIFIER G4double SecondaryThreshold();
  FQUALIFIER void SetSecondaryThreshold(G4double val);
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
  int fThreadId;
  curandState* fDevStates;

  GPElement* fCurrentElement;
  G4int      nsec;
  G4double   xsec[kElementMax];

  G4double energyThresholdLPM;
  G4double klpm,kp;

  G4double particleMass;
  G4double kinEnergy;
  G4double totalEnergy;
  G4double currentZ;
  G4double densityFactor;
  G4double densityCorr;

  G4double lpmEnergy;
  G4double xiLPM;
  G4double phiLPM;
  G4double gLPM;

  G4double z13;
  G4double z23;
  G4double lnZ;
  G4double Fel;
  G4double Finel;
  G4double fCoulomb;
  G4double fMax; 

  G4double lowKinEnergy;

  //secondary from the model
  GXTrack theGamma;

  //SeltzerBergerModel
  GPPhysics2DVector* dataSB;

  //G4VEmModel
  GPVParticleChange* pParticleChange;
  G4double        energyLimitModels;
  G4double        secondaryThreshold;
  G4double        lowLimit;
  G4double        highLimit;
  G4double        eMinActive;
  G4double        eMaxActive;

};

#endif

