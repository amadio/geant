
#include "UserPhysicsList.h"

#include "SystemOfUnits.h"

#include "PhysicsProcess.h"

#include "Particle.h"

#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

#include "ElectronIonizationProcess.h"
#include "MollerBhabhaIonizationModel.h"
#include "ElectronBremsstrahlungProcess.h"
#include "SeltzerBergerBremsModel.h"
#include "RelativisticBremsModel.h"


namespace userapplication {

  UserPhysicsList::UserPhysicsList(const std::string &name) : geantphysics::PhysicsList(name) {}
  UserPhysicsList::~UserPhysicsList() {}

void UserPhysicsList::Initialize() {
  // get the partcile table and loop over it
  std::vector<geantphysics::Particle*> pTable = geantphysics::Particle::GetTheParticleTable();
  for (unsigned int i=0; i<pTable.size(); ++i) {
    geantphysics::Particle *particle = pTable[i];
    if (particle==geantphysics::Electron::Definition()) {
      //std::cout<<"  ELECTRON" <<std::endl;
      //
      // create ionization process for e- with 1 model:
     //
      geantphysics::EMPhysicsProcess *eIoniProc = new geantphysics::ElectronIonizationProcess("eIoni");
      // create the Moller-Bhabha model for ionization i.e. for e- + e- -> e- + e- intercation
      geantphysics::EMModel          *eMBModel  = new geantphysics::MollerBhabhaIonizationModel(true);
      // set min/max energies of the model
      eMBModel->SetLowEnergyUsageLimit (  1.0*geant::keV);
      eMBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add the model to the process
      eIoniProc->AddModel(eMBModel);
      //
      // add the process to the e- particle
      AddProcessToPartcile(particle, eIoniProc);
      //
      // create bremsstrahlung process for e- with 2 models:
      //
      geantphysics::EMPhysicsProcess *eBremProc = new geantphysics::ElectronBremsstrahlungProcess("eBrem");
      // create a SeltzerBergerBremsModel for e-
      geantphysics::EMModel          *eSBModel  = new geantphysics::SeltzerBergerBremsModel(true);
      // set min/max energies of the model
      eSBModel->SetLowEnergyUsageLimit (1.0*geant::keV);
      eSBModel->SetHighEnergyUsageLimit(1.0*geant::GeV);
      // how to inactivate this model in a given region i.e. region with index 1
      // active regions for a model are set based on their process active regions + user requested inactive regions
      //eSBModel->AddToUserRequestedInActiveRegions(1);
      //
      // add this model to the process
      eBremProc->AddModel(eSBModel);
      //
      // create a RelativisticBremsModel for e-
      geantphysics::EMModel          *eRelBModel = new geantphysics::RelativisticBremsModel();
      // set min/max energies of the model
      eRelBModel->SetLowEnergyUsageLimit (  1.0*geant::GeV);
      eRelBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add this model to the process
      eBremProc->AddModel(eRelBModel);
      //
      // add the process to the e- particle
      AddProcessToPartcile(particle, eBremProc);
    }
    if (particle==geantphysics::Positron::Definition()) {
      //std::cout<<"  Positron" <<std::endl;
      //
      // create ionization process for e+ with 1 model:
      //
      geantphysics::EMPhysicsProcess *eIoniProc = new geantphysics::ElectronIonizationProcess("eIoni");
      // create the Moller-Bhabha model for ionization i.e. for e+ + e- -> e+ + e- intercation
      geantphysics::EMModel          *eMBModel  = new geantphysics::MollerBhabhaIonizationModel(false);
      // set min/max energies of the model
      eMBModel->SetLowEnergyUsageLimit (  1.0*geant::keV);
      eMBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add the model to the process
      eIoniProc->AddModel(eMBModel);
      // add the process to the e+ particle
      AddProcessToPartcile(particle, eIoniProc);
      //
      // create bremsstrahlung process for e+ with 2 models:
      //
      geantphysics::EMPhysicsProcess *eBremProc = new geantphysics::ElectronBremsstrahlungProcess("eBrem");
      // create a SeltzerBergerBremsModel for e-
      geantphysics::EMModel          *eSBModel  = new geantphysics::SeltzerBergerBremsModel(false);
      // set min/max energies of the model
      eSBModel->SetLowEnergyUsageLimit (1.0*geant::keV);
      eSBModel->SetHighEnergyUsageLimit(1.0*geant::GeV);
      // how to inactivate this model in a given region i.e. region with index 1
      // active regions for a model are set based on their process active regions + user requested inactive regions
      //eSBModel->AddToUserRequestedInActiveRegions(1);
      //
      // add this model to the process
      eBremProc->AddModel(eSBModel);
      //
      // create a RelativisticBremsModel for e+
      geantphysics::EMModel          *eRelBModel = new geantphysics::RelativisticBremsModel();
      // set min/max energies of the model
      eRelBModel->SetLowEnergyUsageLimit (  1.0*geant::GeV);
      eRelBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add this model to the process
      eBremProc->AddModel(eRelBModel);
      //
      // add the process to the e+ particle
      AddProcessToPartcile(particle, eBremProc);
    }
  }
}


}  // userapplication
