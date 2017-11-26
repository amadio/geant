
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
//#include "G4PhotoElectricEffect.hh"
//#include "G4RayleighScattering.hh"

#include "G4eMultipleScattering.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
//#include "G4eplusAnnihilation.hh"

#include "G4EmParameters.hh"
#include "G4MscStepLimitType.hh"


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
  param->SetMscStepLimitType(fUseSafetyPlus);
  param->SetMscSkin(3);
  param->SetMscRangeFactor(0.1);
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
      ph->RegisterProcess(new G4GammaConversion, particle);
    } else if (particleName =="e-") {
//      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      G4eMultipleScattering* msc         = new G4eMultipleScattering;
      G4GoudsmitSaundersonMscModel* msc1 = new G4GoudsmitSaundersonMscModel();
      msc->AddEmModel(0, msc1);
      ph->RegisterProcess(msc,particle);
      //
      G4eIonisation* eIoni = new G4eIonisation();
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
    } else if (particleName=="e+") {
//      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      G4eMultipleScattering* msc         = new G4eMultipleScattering;
      G4GoudsmitSaundersonMscModel* msc1 = new G4GoudsmitSaundersonMscModel();
      msc->AddEmModel(0, msc1);
      ph->RegisterProcess(msc,particle);
      //
      G4eIonisation* eIoni = new G4eIonisation();
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      //
//      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
    }
  }
}
