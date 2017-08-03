
#include "UserPhysicsList.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

#include "PhysicsProcess.h"

#include "Particle.h"
#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

#include "Proton.h"
#include "Neutron.h"
#include "PionPlus.h"
#include "PionMinus.h"
#include "PionZero.h"
#include "KaonPlus.h"
#include "KaonMinus.h"
#include "KaonZero.h"
#include "KaonShort.h"
#include "KaonLong.h"


#include "ElectronIonizationProcess.h"
#include "MollerBhabhaIonizationModel.h"

#include "ElectronBremsstrahlungProcess.h"
#include "SeltzerBergerBremsModel.h"
#include "RelativisticBremsModel.h"

#include "ComptonScatteringProcess.h"
#include "KleinNishinaComptonModel.h"

#include "GammaConversionProcess.h"
#include "BetheHeitlerPairModel.h"

#include "MSCProcess.h"
#include "MSCModel.h"
#include "GSMSCModel.h"

#include "StepMaxProcess.h"

#include "ElasticScatteringProcess.h"
#include "DiffuseElasticModel.h"
#include "GlauberGribovElasticXsc.h"

namespace userapplication {


UserPhysicsList::UserPhysicsList(const std::string &name) : geantphysics::PhysicsList(name) {
  fMSCSteppingAlgorithm = geantphysics::MSCSteppingAlgorithm::kUseSaftey; // opt0 step limit type
  fStepMaxValue         = geantphysics::PhysicsProcess::GetAVeryLargeValue();
}

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
      geantphysics::EMPhysicsProcess *eIoniProc = new geantphysics::ElectronIonizationProcess("e-Ioni");
      // create the Moller-Bhabha model for ionization i.e. for e- + e- -> e- + e- intercation
      geantphysics::EMModel          *eMBModel  = new geantphysics::MollerBhabhaIonizationModel(true);
      // set min/max energies of the model
      eMBModel->SetLowEnergyUsageLimit (  1.0*geant::keV);
      eMBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
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
      AddProcessToParticle(particle, eBremProc);
      //
      // create MSC process
      geantphysics::EMPhysicsProcess *eMSCProc  = new geantphysics::MSCProcess("e-msc");
      // create GS-msc model, set min/max usage limits
      geantphysics::GSMSCModel       *gsMSCModel = new geantphysics::GSMSCModel();
      gsMSCModel->SetMSCSteppingAlgorithm(fMSCSteppingAlgorithm);
      gsMSCModel->SetLowEnergyUsageLimit(100.*geant::eV);
      gsMSCModel->SetHighEnergyUsageLimit(100.*geant::TeV);
      eMSCProc->AddModel(gsMSCModel);
      // add process to particle
      AddProcessToParticle(particle, eMSCProc);

      //
      // Create and add the special user process
      //
      StepMaxProcess *stepMaxProc = new StepMaxProcess();
      stepMaxProc->SetMaxStep(fStepMaxValue);
      AddProcessToParticle(particle, stepMaxProc);
    }
    if (particle==geantphysics::Positron::Definition()) {
      //std::cout<<"  Positron" <<std::endl;
      //
      // create ionization process for e+ with 1 model:
      //
      geantphysics::EMPhysicsProcess *eIoniProc = new geantphysics::ElectronIonizationProcess("e+Ioni");
      // create the Moller-Bhabha model for ionization i.e. for e+ + e- -> e+ + e- intercation
      geantphysics::EMModel          *eMBModel  = new geantphysics::MollerBhabhaIonizationModel(false);
      // set min/max energies of the model
      eMBModel->SetLowEnergyUsageLimit (  1.0*geant::keV);
      eMBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add the model to the process
      eIoniProc->AddModel(eMBModel);
      // add the process to the e+ particle
      AddProcessToParticle(particle, eIoniProc);
      //
      // create bremsstrahlung process for e+ with 2 models:
      //
      geantphysics::EMPhysicsProcess *eBremProc = new geantphysics::ElectronBremsstrahlungProcess("e+Brem");
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
      AddProcessToParticle(particle, eBremProc);
      //
      // create MSC process
      geantphysics::EMPhysicsProcess *eMSCProc   = new geantphysics::MSCProcess("e+msc");
      // create GS-msc model, set min/max usage limits
      geantphysics::GSMSCModel       *gsMSCModel = new geantphysics::GSMSCModel(false); // for e+
      gsMSCModel->SetMSCSteppingAlgorithm(fMSCSteppingAlgorithm);
      gsMSCModel->SetLowEnergyUsageLimit(100.*geant::eV);
      gsMSCModel->SetHighEnergyUsageLimit(100.*geant::TeV);
      gsMSCModel->SetOptionPWAScreening(true);
      eMSCProc->AddModel(gsMSCModel);
      // add process to particle
      AddProcessToParticle(particle, eMSCProc);

      //
      // Create and add the special user process
      //
      StepMaxProcess *stepMaxProc = new StepMaxProcess();
      stepMaxProc->SetMaxStep(fStepMaxValue);
      AddProcessToParticle(particle, stepMaxProc);
    }
    if (particle==geantphysics::Gamma::Definition()) {
      // create compton scattering process for gamma with 1 model:
      //
      geantphysics::EMPhysicsProcess *comptProc = new geantphysics::ComptonScatteringProcess();
      // create the Klein-Nishina model for Compton scattering i.e. for g + e- -> g + e- intercation
      geantphysics::EMModel          *kncModel  = new geantphysics::KleinNishinaComptonModel();
      // set min/max energies of the model
      kncModel->SetLowEnergyUsageLimit (100.0*geant::eV);
      kncModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add the model to the process
      comptProc->AddModel(kncModel);
      //
      // add the process to the gamma particle
      AddProcessToParticle(particle, comptProc);
      //
      // create gamma conversion process for gamma with 1 model:
      //
      geantphysics::EMPhysicsProcess *convProc = new geantphysics::GammaConversionProcess();
      // create the Bethe-Heitler model for pair production i.e. for g + A -> e- + e+ intercation
      geantphysics::EMModel           *bhModel = new geantphysics::BetheHeitlerPairModel();
      // set min/max energies of the model
      bhModel->SetLowEnergyUsageLimit (  2.0*geant::kElectronMassC2);
      // the parametrized cross sections works only up t0 80-90 GeV but we will use it now up to 1 TeV
      // it will be changed when we will have the high-energy model
      bhModel->SetHighEnergyUsageLimit(  1.0*geant::TeV);
      // add the model to the process
      convProc->AddModel(bhModel);
      //
      // add the process to the gamma particle
      AddProcessToParticle(particle, convProc);
    }
    if (particle==geantphysics::Proton::Definition() || 
	particle==geantphysics::Neutron::Definition() ||
	particle==geantphysics::PionPlus::Definition() ||
	particle==geantphysics::PionMinus::Definition() ||
	particle==geantphysics::PionZero::Definition() ||
	particle==geantphysics::KaonPlus::Definition() ||
	particle==geantphysics::KaonMinus::Definition() ||
	particle==geantphysics::KaonZero::Definition() ||
	particle==geantphysics::KaonShort::Definition() ||
	particle==geantphysics::KaonLong::Definition()) {
      // create hadronic elastic process for proton:
      //
      geantphysics::HadronicProcess *helProc = new geantphysics::ElasticScatteringProcess();
      // create the diffuse elastic model for elastic scattering
      geantphysics::HadronicFinalStateModel *diffelModel  = new geantphysics::DiffuseElasticModel();
      // create the cross sections
      geantphysics::HadronicCrossSection *ggElasticXS = new geantphysics::GlauberGribovElasticXsc();
      
      // set min/max energies of the model
      diffelModel->SetLowEnergyUsageLimit (100.0*geant::eV);
      diffelModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
      // add the model to the process
      helProc->AddModel(diffelModel);
      // add the cross-sections to the process
      helProc->AddCrossSection(ggElasticXS);
      //
      // add the process to the gamma particle
      AddProcessToParticle(particle, helProc);
    }
  }
}

void  UserPhysicsList::SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm stepping) {
  fMSCSteppingAlgorithm = stepping;
}

void  UserPhysicsList::SetStepMaxValue(double val) { fStepMaxValue = val; }


}  // userapplication
