
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

#include "ElasticScatteringProcess.h"
#include "DiffuseElasticModel.h"
#include "GlauberGribovElasticXsc.h"

namespace userapplication {

  UserPhysicsList::UserPhysicsList(const std::string &name) : geantphysics::PhysicsList(name) {}
  UserPhysicsList::~UserPhysicsList() {}

  void UserPhysicsList::Initialize() {
    // get the partcile table and loop over it
    std::vector<geantphysics::Particle*> pTable = geantphysics::Particle::GetTheParticleTable();
    for (unsigned int i=0; i<pTable.size(); ++i) {
      geantphysics::Particle *particle = pTable[i];
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


}  // userapplication
