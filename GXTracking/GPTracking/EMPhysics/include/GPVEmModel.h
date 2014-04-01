#ifndef GPVEmModel_h
#define GPVEmModel_h 1

#include "GPTypeDef.h"
#include "GPMaterial.h"
#include "GPElement.h"

#include "GPPhysicsTable.h"
#include "GPConstants.h"
#include "GPVParticleChange.h"
#include "GXTrack.h"

class GPVEmModel
{
public:

  FQUALIFIER GPVEmModel(curandState* devStates, int threadId);
  FQUALIFIER ~GPVEmModel();

  FQUALIFIER void SampleSecondaries(GXTrack* track,
				    GPMaterial* material,
				    G4double tmin = 0.0,
				    G4double tmax = DBL_MAX);

  FQUALIFIER GXTrack& GetSecondary();

  // main method to compute cross section per Volume
  FQUALIFIER G4double G4VEmModel_CrossSectionPerVolume(GPMaterial*,
						       G4double kineticEnergy,
						       G4double cutEnergy = 0.0,
						       G4double maxEnergy = DBL_MAX);

  //@@@ main method to compute cross section per Volume of a derived Class
  FQUALIFIER G4double CrossSectionPerVolume(GPMaterial*,
					    G4double kineticEnergy,
					    G4double cutEnergy = 0.0,
					    G4double maxEnergy = DBL_MAX);


  // to select atom cross section per volume is recomputed for each element 
  FQUALIFIER GPElement* SelectRandomAtom(GPMaterial* material,
					 G4double kineticEnergy,
					 G4double cutEnergy = 0.0,
					 G4double maxEnergy = DBL_MAX);
  
  // main method to compute cross section per atom
  FQUALIFIER G4double ComputeCrossSectionPerAtom(G4double kinEnergy, 
						 G4double Z, 
						 G4double A = 0., /* amu */
						 G4double cutEnergy = 0.0,
						 G4double maxEnergy = DBL_MAX);
				      				     
  // main method to compute dEdx
  FQUALIFIER G4double ComputeDEDXPerVolume(GPMaterial*,
  					   G4double kineticEnergy,
  					   G4double cutEnergy = DBL_MAX);
  
  // Compute effective ion charge square
  FQUALIFIER G4double ChargeSquareRatio(GXTrack* track);

  // Compute effective ion charge square
  FQUALIFIER G4double GetChargeSquareRatio(G4double q);

  // Compute ion charge 
  FQUALIFIER G4double GetParticleCharge();

  // value which may be tabulated (by default cross section)
  FQUALIFIER G4double Value(GPMaterial* material,
			    G4double kineticEnergy);

  // threshold for zero value 
  FQUALIFIER G4double MinPrimaryEnergy();

  // kinematically allowed max kinetic energy of a secondary
  FQUALIFIER G4double MaxSecondaryEnergy(G4double kineticEnergy);  

  // initilisation at run time for a given material
  FQUALIFIER void SetupForMaterial(GPMaterial*,
				   G4double kineticEnergy);

  FQUALIFIER void SetParticleChange(GPVParticleChange* p);

  FQUALIFIER void SetCrossSectionTable(GPPhysicsTable* p);

  FQUALIFIER GPElement* GetCurrentElement();

  FQUALIFIER void SetCurrentElement(GPElement* elm);

  FQUALIFIER G4double MaxSecondaryKinEnergy(G4double kineticEnergy);

  // dEdx per unit length
  FQUALIFIER G4double ComputeDEDX(GPMaterial* material,
				  G4double kineticEnergy,
				  G4double cutEnergy = DBL_MAX);

  // cross section per volume
  FQUALIFIER G4double CrossSection(GPMaterial* material,
				   G4double kineticEnergy,
				   G4double cutEnergy = 0.0,
				   G4double maxEnergy = DBL_MAX);

  // compute mean free path via cross section per volume
  FQUALIFIER G4double ComputeMeanFreePath(G4double kineticEnergy,
					  GPMaterial* material,    
					  G4double cutEnergy = 0.0,
					  G4double maxEnergy = DBL_MAX);

  // generic cross section per element
  FQUALIFIER G4double ComputeCrossSectionPerAtom(GPElement* elm,
						 G4double kinEnergy, 
						 G4double cutEnergy = 0.0,
						 G4double maxEnergy = DBL_MAX);
  
  //  FQUALIFIER G4VEmAngularDistribution* GetAngularDistribution();

  //  FQUALIFIER void SetAngularDistribution(G4VEmAngularDistribution*);

  FQUALIFIER G4double HighEnergyLimit();

  FQUALIFIER G4double LowEnergyLimit();

  FQUALIFIER G4double HighEnergyActivationLimit();

  FQUALIFIER G4double LowEnergyActivationLimit();

  FQUALIFIER G4double PolarAngleLimit();

  FQUALIFIER G4double SecondaryThreshold();

  FQUALIFIER G4bool LPMFlag();

  FQUALIFIER void SetHighEnergyLimit(G4double val);

  FQUALIFIER void SetLowEnergyLimit(G4double val);

  FQUALIFIER void SetActivationHighEnergyLimit(G4double val);

  FQUALIFIER void SetActivationLowEnergyLimit(G4double val);

  FQUALIFIER G4bool IsActive(G4double kinEnergy);

  FQUALIFIER void SetPolarAngleLimit(G4double val);

  FQUALIFIER void SetSecondaryThreshold(G4double val);

  FQUALIFIER void SetLPMFlag(G4bool val);

  FQUALIFIER GPPhysicsTable* GetCrossSectionTable();

  // initialisation of the ParticleChange for the model
  FQUALIFIER GPVParticleChange* GetParticleChange();

  // initialisation of the ParticleChange for the model
  //  G4ParticleChangeForGamma* GetParticleChangeForGamma();

private:

  //curand                  
  int fThreadId;
  curandState* fDevStates;

  G4double        lowLimit;
  G4double        highLimit;
  G4double        eMinActive;
  G4double        eMaxActive;
  G4double        polarAngleLimit;
  G4double        secondaryThreshold;
  G4bool          theLPMflag;

  GPVParticleChange*  pParticleChange;
  GPPhysicsTable*     xSectionTable;
  GPElement*          fCurrentElement;

  G4double  xsec[nElements];

  //the secondary electron
  GXTrack theSecondaryElectron;

};

#endif

