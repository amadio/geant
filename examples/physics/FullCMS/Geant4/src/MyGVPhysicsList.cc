
#include "MyGVPhysicsList.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
//#include "G4RayleighScattering.hh"

#include "G4eMultipleScattering.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4EmParameters.hh"
#include "G4MscStepLimitType.hh"

#include "G4EmConfigurator.hh"
#include "G4LossTableManager.hh"

MyGVPhysicsList::MyGVPhysicsList() : G4VUserPhysicsList() {
  SetDefaultCutValue(1.0);
  SetVerboseLevel(0);
}


MyGVPhysicsList::~MyGVPhysicsList() {}


void MyGVPhysicsList::ConstructParticle() {
   G4Electron::ElectronDefinition();
   G4Positron::PositronDefinition();
   G4Gamma::GammaDefinition();
}


void MyGVPhysicsList::ConstructProcess() {
  // Transportation
  AddTransportation();
  // EM physics
  BuildEMPhysics();
}


void MyGVPhysicsList::BuildEMPhysics() {
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(1);
  // inactivate energy loss fluctuations
  param->SetLossFluctuations(false);
  // inactivate to use cuts as final range
  param->SetUseCutAsFinalRange(false);
  //
  // MSC options:
  param->SetMscStepLimitType(fUseSafety);
  param->SetMscSkin(3);
  param->SetMscRangeFactor(0.2);  // default EM-opt1 value
  G4LossTableManager::Instance();
  //
  // Add standard EM physics processes to e-/e+ and gamma that GeantV has
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  auto aParticleIterator  = GetParticleIterator();
  aParticleIterator->reset();
  while((*aParticleIterator)()) {
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String particleName          = particle->GetParticleName();
    if (particleName=="gamma") {
//      ph->RegisterProcess(new G4PhotoElectricEffect, particle);
      ph->RegisterProcess(new G4ComptonScattering(), particle);
      ph->RegisterProcess(new G4GammaConversion(), particle);
      G4double LivermoreLowEnergyLimit  = 1.*eV;
      G4double LivermoreHighEnergyLimit = 1.*TeV;
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = new G4LivermorePhotoElectricModel();
      theLivermorePhotoElectricModel->SetLowEnergyLimit(LivermoreLowEnergyLimit);
      theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
      ph->RegisterProcess(thePhotoElectricEffect, particle);
    } else if (particleName =="e-") {
//      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      G4eMultipleScattering* msc = new G4eMultipleScattering;
      // 2 msc models:
      // = 1. with range-factor = 0.2  (EM-opt1) settings and active everywhere but in HcalRegion (#11)
      // = 2. with range-factor = 0.06 (EM-opt0) settings and active only in HcalRegion (#11)
      G4GoudsmitSaundersonMscModel* mscOpt1 = new G4GoudsmitSaundersonMscModel();
      G4GoudsmitSaundersonMscModel* mscOpt0 = new G4GoudsmitSaundersonMscModel();
      // set the 2. msc model to be opt0 and lock it to prevet update at init.
      mscOpt0->SetRangeFactor(0.06); // all others are the default
      mscOpt0->SetSkin(3);           // not important just for consistency
      mscOpt0->SetLocked(true);
      // add the default opt1 model to the process and register the process for the particle
      msc->AddEmModel(0, mscOpt1);
      ph->RegisterProcess(msc,particle);
      // Add the extra mscOpt0 model to e- and msc-process used in the "HcalRegion"
      G4EmConfigurator* emConfig = G4LossTableManager::Instance()->EmConfigurator();
      emConfig->SetExtraEmModel("e-", "msc", mscOpt0, "HcalRegion");
      //
      //
      G4eIonisation* eIoni = new G4eIonisation();
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
    } else if (particleName=="e+") {
//      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      G4eMultipleScattering* msc = new G4eMultipleScattering;
      // 2 msc models:
      // = 1. with range-factor = 0.2  (EM-opt1) settings and active everywhere but in HcalRegion (#11)
      // = 2. with range-factor = 0.06 (EM-opt0) settings and active only in HcalRegion (#11)
      G4GoudsmitSaundersonMscModel* mscOpt1 = new G4GoudsmitSaundersonMscModel();
      G4GoudsmitSaundersonMscModel* mscOpt0 = new G4GoudsmitSaundersonMscModel();
      // set the 2. msc model to be opt0 and lock it to prevet update at init.
      mscOpt0->SetRangeFactor(0.06); // all others are the default
      mscOpt0->SetSkin(3);           // not important just for consistency
      mscOpt0->SetLocked(true);
      // add the default opt1 model to the process and register the process for the particle
      msc->AddEmModel(0, mscOpt1);
      ph->RegisterProcess(msc,particle);
      // Add the extra mscOpt0 model to e+ and msc-process used in the "HcalRegion"
      G4EmConfigurator* emConfig = G4LossTableManager::Instance()->EmConfigurator();
      emConfig->SetExtraEmModel("e+", "msc", mscOpt0, "HcalRegion");
      //
      //
      G4eIonisation* eIoni = new G4eIonisation();
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      //
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
    }
  }
}
