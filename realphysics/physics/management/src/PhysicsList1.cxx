
#include "PhysicsList1.h"

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

#include "ComptonScatteringProcess.h"
#include "KleinNishinaComptonModel.h"


namespace geantphysics {

  PhysicsList1::PhysicsList1(const std::string &name) : PhysicsList(name) {}
  PhysicsList1::~PhysicsList1() {}

void PhysicsList1::Initialize() {
  // get the partcile table and loop over it
  std::vector<Particle*> pTable = Particle::GetTheParticleTable();
  for (unsigned int i=0; i<pTable.size(); ++i) {
    Particle *particle = pTable[i];
    if (particle==Electron::Definition()) {
      //std::cout<<"  ELECTRON" <<std::endl;
      //
      // create ionization process for e- with 1 model:
     //
      EMPhysicsProcess *eIoniProc = new ElectronIonizationProcess("eIoni");
      // create the Moller-Bhabha model for ionization i.e. for e- + e- -> e- + e- intercation
      EMModel          *eMBModel  = new MollerBhabhaIonizationModel(true);
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
      EMPhysicsProcess *eBremProc = new ElectronBremsstrahlungProcess("eBrem");
      // create a SeltzerBergerBremsModel for e-
      EMModel          *eSBModel  = new SeltzerBergerBremsModel(true);
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
      EMModel          *eRelBModel = new RelativisticBremsModel();
      // set min/max energies of the model
      eRelBModel->SetLowEnergyUsageLimit (  1.0*geant::GeV);
      eRelBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add this model to the process
      eBremProc->AddModel(eRelBModel);
      //
      // add the process to the e- particle
      AddProcessToPartcile(particle, eBremProc);
    }
    if (particle==Positron::Definition()) {
      //std::cout<<"  Positron" <<std::endl;
      //
      // create ionization process for e+ with 1 model:
      //
      EMPhysicsProcess *eIoniProc = new ElectronIonizationProcess("eIoni");
      // create the Moller-Bhabha model for ionization i.e. for e+ + e- -> e+ + e- intercation
      EMModel          *eMBModel  = new MollerBhabhaIonizationModel(false);
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
      EMPhysicsProcess *eBremProc = new ElectronBremsstrahlungProcess("eBrem");
      // create a SeltzerBergerBremsModel for e-
      EMModel          *eSBModel  = new SeltzerBergerBremsModel(false);
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
      EMModel          *eRelBModel = new RelativisticBremsModel();
      // set min/max energies of the model
      eRelBModel->SetLowEnergyUsageLimit (  1.0*geant::GeV);
      eRelBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add this model to the process
      eBremProc->AddModel(eRelBModel);
      //
      // add the process to the e+ particle
      AddProcessToPartcile(particle, eBremProc);
    }
    if (particle==Gamma::Definition()) {
      // create compton scattering process for e- with 1 model:
      //
      EMPhysicsProcess *comptProc = new ComptonScatteringProcess();
      // create the Klein-Nishina model for Compton scattering i.e. for g + e- -> g + e- intercation
      EMModel          *kncModel = new KleinNishinaComptonModel();
      // set min/max energies of the model
      kncModel->SetLowEnergyUsageLimit (100.0*geant::eV);
      kncModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add the model to the process
      comptProc->AddModel(kncModel);
      //
      // add the process to the gamma particle
      AddProcessToPartcile(particle, comptProc);
    }

  }
}


}  // namespace geantphysics
