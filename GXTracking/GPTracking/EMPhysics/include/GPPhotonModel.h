#ifndef GPPhotonModel_h
#define GPPhotonModel_h 1

#include "GPTypeDef.h"
#include "GPEmProcessType.h"
#include "GPMaterial.h"
#include "GPElement.h"
#include "GPSandiaTable.h"

#include "GPPhysicsTable.h"
#include "GPConstants.h"
#include "GPVParticleChange.h"
#include "GXTrack.h"

#include "GPRandom.h"

class GPPhotonModel
{
public:

  //---------------------------------------------------------------------------
  // G4VEmModel 
  //---------------------------------------------------------------------------

  FQUALIFIER GPPhotonModel(curandState* devStates, int threadId);

  FQUALIFIER ~GPPhotonModel();

  FQUALIFIER void SetProcessType(GPEmProcessType aType);

  FQUALIFIER void SampleSecondaries(GXTrack* track,
				    GPMaterial* material,
				    G4double tmin = 0.0,
				    G4double tmax = DBL_MAX);

  // main method to compute cross section per Volume
  FQUALIFIER 
  G4double G4VEmModel_CrossSectionPerVolume(GPMaterial* material,
					    G4double kineticEnergy,
					    G4double cutEnergy = 0.0,
					    G4double maxEnergy = DBL_MAX);

  //@@@ main method to compute cross section per Volume of a derived Class
  FQUALIFIER G4double CrossSectionPerVolume(GPMaterial* material,
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
  FQUALIFIER G4double ComputeDEDXPerVolume(GPMaterial* material,
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

  //---------------------------------------------------------------------------
  // G4KleinNishinaCompton for G4ComptonScattering
  //---------------------------------------------------------------------------

  FQUALIFIER G4double G4KleinNishinaCompton_ComputeCrossSectionPerAtom(
                                                 G4double kinEnergy, 
						 G4double Z);
  FQUALIFIER 
  void G4KleinNishinaCompton_SampleSecondaries(GXTrack* aDynamicGamma);

  //---------------------------------------------------------------------------
  // G4PEEffectFluoModel for G4PhotoElectricEffect
  //---------------------------------------------------------------------------

  FQUALIFIER
  G4double G4PEEffectFluoModel_ComputeCrossSectionPerAtom(G4double energy,
							  G4double Z);

  FQUALIFIER
  G4double G4PEEffectFluoModel_CrossSectionPerVolume(GPMaterial* material,
						     G4double energy);
  FQUALIFIER 
  void G4PEEffectFluoModel_SampleSecondaries(GPMaterial* aMaterial,
                                             GXTrack* aDynamicPhoton);

  FQUALIFIER 
  GPThreeVector SauterGavrilaAngularDistribution_SampleDirection(GXTrack* 
								 aDynamicGamma);

  FQUALIFIER
  void SetaSandiaTable(GPSandiaTable* table);

  //---------------------------------------------------------------------------
  // G4BetheHeitlerModel for G4GammaConversion
  //---------------------------------------------------------------------------

  FQUALIFIER
  void G4BetheHeitlerModel_SampleSecondaries(GPMaterial* aMaterial,
					     GXTrack* aDynamicGamma);

  FQUALIFIER
  G4double G4BetheHeitlerModel_ScreenFunction1(G4double ScreenVariable);

  FQUALIFIER 
  G4double G4BetheHeitlerModel_ScreenFunction2(G4double ScreenVariable);

  FQUALIFIER
  G4double G4BetheHeitlerModel_ComputeCrossSectionPerAtom(G4double GammaEnergy, 
							  G4double Z);

  //-------------------------------------------------------------------------
  // Utility Functions for Seoncdaries
  //-------------------------------------------------------------------------
  FQUALIFIER GXTrack& GetSecondaryElectron();

  FQUALIFIER GXTrack& GetSecondaryPositron();

  FQUALIFIER void FillSecondary(GXTrack* secondary, GXTrack* track, 
				GPThreeVector direction,
				G4double energy, G4double charge);

private:

  //curand                  
  int fThreadId;
  curandState* fDevStates;

  //Process
  GPEmProcessType              theProcessType;

  //GPVEmModel
  G4double        lowLimit;
  G4double        highLimit;
  G4double        eMinActive;
  G4double        eMaxActive;
  G4double        polarAngleLimit;
  G4double        secondaryThreshold;
  G4bool          theLPMflag;

  GPVParticleChange*  fParticleChange;
  GPPhysicsTable*     xSectionTable;
  GPElement*          fCurrentElement;

  G4double  xsec[nElements];

  //the secondary electron
  GXTrack theSecondaryElectron;
  GXTrack theSecondaryPositron;

  //used in G4PEEffectFluoModel - should move to GPMaterial
  GPSandiaTable* fSandiaTable;

};

#endif

