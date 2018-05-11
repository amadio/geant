#include "Geant/PositronTo2GammaModel.h"
#include "Geant/GSMSCModelSimplified.h"
#include "TestEm3PhysicsList.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

#include "Geant/PhysicsProcess.h"

#include "Geant/Particle.h"
#include "Geant/Electron.h"
#include "Geant/Positron.h"
#include "Geant/Gamma.h"

#include "Geant/ElectronIonizationProcess.h"
#include "Geant/MollerBhabhaIonizationModel.h"

#include "Geant/ElectronBremsstrahlungProcess.h"
#include "Geant/SeltzerBergerBremsModel.h"
#include "Geant/RelativisticBremsModel.h"

#include "Geant/PositronAnnihilationProcess.h"

#include "Geant/ComptonScatteringProcess.h"
#include "Geant/KleinNishinaComptonModel.h"

#include "Geant/GammaConversionProcess.h"
#include "Geant/BetheHeitlerPairModel.h"
#include "Geant/RelativisticPairModel.h"

#include "Geant/GammaPhotoElectricProcess.h"
#include "Geant/SauterGavrilaPhotoElectricModel.h"

#include "Geant/MSCProcess.h"
#include "Geant/MSCModel.h"
#include "Geant/GSMSCModel.h"

#include "TestEm3FastSimProcess.h"

namespace userapplication {

TestEm3PhysicsList::TestEm3PhysicsList(const std::string &name, bool vector)
    : geantphysics::PhysicsList(name), fVectorized(vector)
{
  fMSCSteppingAlgorithm = geantphysics::MSCSteppingAlgorithm::kUseSaftey; // opt0 step limit type
}

TestEm3PhysicsList::~TestEm3PhysicsList() {}

void TestEm3PhysicsList::Initialize()
{
  // get the partcile table and loop over it
  std::vector<geantphysics::Particle *> pTable = geantphysics::Particle::GetTheParticleTable();
  for (unsigned int i = 0; i < pTable.size(); ++i) {
    geantphysics::Particle *particle = pTable[i];
    if (particle == geantphysics::Electron::Definition()) {
      // std::cout<<"  ELECTRON" <<std::endl;
      //
      // create ionization process for e- with 1 model:
      //
      geantphysics::EMPhysicsProcess *eIoniProc = new geantphysics::ElectronIonizationProcess("e-Ioni");
      // create the Moller-Bhabha model for ionization i.e. for e- + e- -> e- + e- intercation
      geantphysics::EMModel *eMBModel = new geantphysics::MollerBhabhaIonizationModel(true);
      eMBModel->SetBasketizable(fVectorized);
      // set min/max energies of the model
      eMBModel->SetLowEnergyUsageLimit(1.0 * geant::units::keV);
      eMBModel->SetHighEnergyUsageLimit(100.0 * geant::units::TeV);
      // add the model to the process
      eIoniProc->AddModel(eMBModel);
      //
      // add the process to the e- particle
      AddProcessToParticle(particle, eIoniProc);
      //
      // create bremsstrahlung process for e- with 2 models:
      //
      geantphysics::EMPhysicsProcess *eBremProc = new geantphysics::ElectronBremsstrahlungProcess("e-Brem");
      // create a SeltzerBergerBremsModel for e-
      geantphysics::EMModel *eSBModel = new geantphysics::SeltzerBergerBremsModel(true);
      eSBModel->SetBasketizable(fVectorized);
      // set min/max energies of the model
      eSBModel->SetLowEnergyUsageLimit(1.0 * geant::units::keV);
      eSBModel->SetHighEnergyUsageLimit(1.0 * geant::units::GeV);
      // how to inactivate this model in a given region i.e. region with index 1
      // active regions for a model are set based on their process active regions + user requested inactive regions
      // eSBModel->AddToUserRequestedInActiveRegions(1);
      //
      // add this model to the process
      eBremProc->AddModel(eSBModel);
      //
      // create a RelativisticBremsModel for e-
      geantphysics::EMModel *eRelBModel = new geantphysics::RelativisticBremsModel();
      eRelBModel->SetBasketizable(fVectorized);
      // set min/max energies of the model
      eRelBModel->SetLowEnergyUsageLimit(1.0 * geant::units::GeV);
      eRelBModel->SetHighEnergyUsageLimit(100.0 * geant::units::TeV);
      // add this model to the process
      eBremProc->AddModel(eRelBModel);
      //
      // add the process to the e- particle
      AddProcessToParticle(particle, eBremProc);
      //
      // create MSC process
      geantphysics::EMPhysicsProcess *eMSCProc = new geantphysics::MSCProcess("e-msc");
      // create GS-msc model, set min/max usage limits
      geantphysics::MSCModel *gsMSCModel;
      if (fVectorized) {
        gsMSCModel = new geantphysics::GSMSCModelSimplified(true);
      } else {
        gsMSCModel = new geantphysics::GSMSCModel(true);
      }
      gsMSCModel->SetRangeFactor(0.06);
      gsMSCModel->SetMSCSteppingAlgorithm(fMSCSteppingAlgorithm);
      gsMSCModel->SetLowEnergyUsageLimit(100. * geant::units::eV);
      gsMSCModel->SetHighEnergyUsageLimit(100. * geant::units::TeV);
      eMSCProc->AddModel(gsMSCModel);
      // add process to particle
      AddProcessToParticle(particle, eMSCProc);
    }
    if (particle == geantphysics::Positron::Definition()) {
      // std::cout<<"  Positron" <<std::endl;
      //
      // create ionization process for e+ with 1 model:
      //
      geantphysics::EMPhysicsProcess *eIoniProc = new geantphysics::ElectronIonizationProcess("e+Ioni");
      // create the Moller-Bhabha model for ionization i.e. for e+ + e- -> e+ + e- intercation
      geantphysics::EMModel *eMBModel = new geantphysics::MollerBhabhaIonizationModel(false);
      eMBModel->SetBasketizable(fVectorized);
      // set min/max energies of the model
      eMBModel->SetLowEnergyUsageLimit(1.0 * geant::units::keV);
      eMBModel->SetHighEnergyUsageLimit(100.0 * geant::units::TeV);
      // add the model to the process
      eIoniProc->AddModel(eMBModel);
      // add the process to the e+ particle
      AddProcessToParticle(particle, eIoniProc);
      //
      // create bremsstrahlung process for e+ with 2 models:
      //
      geantphysics::EMPhysicsProcess *eBremProc = new geantphysics::ElectronBremsstrahlungProcess("e+Brem");
      // create a SeltzerBergerBremsModel for e-
      geantphysics::EMModel *eSBModel = new geantphysics::SeltzerBergerBremsModel(false);
      eSBModel->SetBasketizable(fVectorized);
      // set min/max energies of the model
      eSBModel->SetLowEnergyUsageLimit(1.0 * geant::units::keV);
      eSBModel->SetHighEnergyUsageLimit(1.0 * geant::units::GeV);
      // how to inactivate this model in a given region i.e. region with index 1
      // active regions for a model are set based on their process active regions + user requested inactive regions
      // eSBModel->AddToUserRequestedInActiveRegions(1);
      //
      // add this model to the process
      eBremProc->AddModel(eSBModel);
      //
      // create a RelativisticBremsModel for e+
      geantphysics::EMModel *eRelBModel = new geantphysics::RelativisticBremsModel();
      eRelBModel->SetBasketizable(fVectorized);
      // set min/max energies of the model
      eRelBModel->SetLowEnergyUsageLimit(1.0 * geant::units::GeV);
      eRelBModel->SetHighEnergyUsageLimit(100.0 * geant::units::TeV);
      // add this model to the process
      eBremProc->AddModel(eRelBModel);
      //
      // add the process to the e+ particle
      AddProcessToParticle(particle, eBremProc);
      //
      // create MSC process
      geantphysics::EMPhysicsProcess *eMSCProc = new geantphysics::MSCProcess("e+msc");
      // create GS-msc model, set min/max usage limits
      geantphysics::MSCModel *gsMSCModel;
      if (fVectorized) {
        gsMSCModel = new geantphysics::GSMSCModelSimplified(false);
      } else {
        gsMSCModel = new geantphysics::GSMSCModel(false);
      }
      gsMSCModel->SetRangeFactor(0.06);
      gsMSCModel->SetMSCSteppingAlgorithm(fMSCSteppingAlgorithm);
      gsMSCModel->SetLowEnergyUsageLimit(100. * geant::units::eV);
      gsMSCModel->SetHighEnergyUsageLimit(100. * geant::units::TeV);
      eMSCProc->AddModel(gsMSCModel);
      // add process to particle
      AddProcessToParticle(particle, eMSCProc);
      //
      // create e+ electron annihilation into 2 gamma process
      geantphysics::PositronAnnihilationProcess *pAnhProc = new geantphysics::PositronAnnihilationProcess();
      geantphysics::PositronTo2GammaModel *pos2Gamma      = new geantphysics::PositronTo2GammaModel();
      pos2Gamma->SetBasketizable(fVectorized);
      pos2Gamma->SetLowEnergyUsageLimit(100. * geant::units::eV);
      pos2Gamma->SetHighEnergyUsageLimit(100. * geant::units::TeV);
      pAnhProc->AddModel(pos2Gamma);
      AddProcessToParticle(particle, pAnhProc);
    }
    if (particle == geantphysics::Gamma::Definition()) {
      // create compton scattering process for gamma with 1 model:
      //
      geantphysics::EMPhysicsProcess *comptProc = new geantphysics::ComptonScatteringProcess();
      // create the Klein-Nishina model for Compton scattering i.e. for g + e- -> g + e- intercation
      geantphysics::EMModel *kncModel = new geantphysics::KleinNishinaComptonModel();
      kncModel->SetBasketizable(fVectorized);
      // set min/max energies of the model
      kncModel->SetLowEnergyUsageLimit(100.0 * geant::units::eV);
      kncModel->SetHighEnergyUsageLimit(100.0 * geant::units::TeV);
      // add the model to the process
      comptProc->AddModel(kncModel);
      //
      // add the process to the gamma particle
      AddProcessToParticle(particle, comptProc);
      //
      // create gamma conversion process for gamma with 1 model:
      //
      geantphysics::EMPhysicsProcess *convProc = new geantphysics::GammaConversionProcess();
      // create the Bethe-Heitler model for pair production i.e. for g + A -> e- + e+ interaction
      geantphysics::EMModel *bhModel = new geantphysics::BetheHeitlerPairModel();
      bhModel->SetBasketizable(fVectorized);
      // set min/max energies of the model
      bhModel->SetLowEnergyUsageLimit(2.0 * geant::units::kElectronMassC2);
      bhModel->SetHighEnergyUsageLimit(80.0 * geant::units::GeV);
      // add the model to the process
      convProc->AddModel(bhModel);
      //
      // create the relativistic model(with LPM) for pair production i.e. for g + A -> e- + e+ interaction
      geantphysics::EMModel *relModel = new geantphysics::RelativisticPairModel();
      relModel->SetBasketizable(fVectorized);
      // set min/max energies of the model
      relModel->SetLowEnergyUsageLimit(80.0 * geant::units::GeV);
      relModel->SetHighEnergyUsageLimit(100.0 * geant::units::TeV);
      // add the model to the process
      convProc->AddModel(relModel);
      //
      // add the process to the gamma particle
      AddProcessToParticle(particle, convProc);
      //
      // create photoelectric effect process for gamma with 1 model:
      //
      geantphysics::EMPhysicsProcess *photoelectricProc = new geantphysics::GammaPhotoElectricProcess();
      // create the Sauter-Gavrila model for photoelectric effect
      geantphysics::EMModel *sgModel = new geantphysics::SauterGavrilaPhotoElectricModel();
      // set min/max energies of the model
      sgModel->SetLowEnergyUsageLimit(1.0 * geant::units::eV);
      sgModel->SetHighEnergyUsageLimit(1.0 * geant::units::TeV);
      // add the model to the process
      photoelectricProc->AddModel(sgModel);
      //
      // add the process to the gamma particle
      AddProcessToParticle(particle, photoelectricProc);

      // addding Fast Sim process for gamma
      geantphysics::FastSimProcess *fastsimproc = new TestEm3FastSimProcess();
      AddProcessToParticle(particle, fastsimproc);
    }
  }
}

void TestEm3PhysicsList::SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm stepping)
{
  fMSCSteppingAlgorithm = stepping;
}

} // namespace userapplication
