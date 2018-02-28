
#include "UserPhysicsList.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

#include "Geant/PhysicsProcess.h"

#include "Geant/Particle.h"
#include "Geant/Electron.h"
#include "Geant/Positron.h"
#include "Geant/Gamma.h"

#include "Geant/Proton.h"
#include "Geant/Neutron.h"
#include "Geant/PionPlus.h"
#include "Geant/PionMinus.h"
#include "Geant/PionZero.h"
#include "Geant/KaonPlus.h"
#include "Geant/KaonMinus.h"
#include "Geant/KaonZero.h"
#include "Geant/KaonShort.h"
#include "Geant/KaonLong.h"

#include "Geant/ElectronIonizationProcess.h"
#include "Geant/MollerBhabhaIonizationModel.h"

#include "Geant/ElectronBremsstrahlungProcess.h"
#include "Geant/SeltzerBergerBremsModel.h"
#include "Geant/RelativisticBremsModel.h"

#include "Geant/ComptonScatteringProcess.h"
#include "Geant/KleinNishinaComptonModel.h"

#include "Geant/GammaConversionProcess.h"
#include "Geant/BetheHeitlerPairModel.h"

#include "Geant/ElasticScatteringProcess.h"
#include "Geant/DiffuseElasticModel.h"
#include "Geant/GlauberGribovElasticXsc.h"

namespace userapplication {

UserPhysicsList::UserPhysicsList(const std::string &name) : geantphysics::PhysicsList(name)
{
}
UserPhysicsList::~UserPhysicsList()
{
}

void UserPhysicsList::Initialize()
{
  // get the partcile table and loop over it
  std::vector<geantphysics::Particle *> pTable = geantphysics::Particle::GetTheParticleTable();
  for (unsigned int i = 0; i < pTable.size(); ++i) {
    geantphysics::Particle *particle = pTable[i];
    if (particle == geantphysics::Proton::Definition() || particle == geantphysics::Neutron::Definition() ||
        particle == geantphysics::PionPlus::Definition() || particle == geantphysics::PionMinus::Definition() ||
        particle == geantphysics::PionZero::Definition() || particle == geantphysics::KaonPlus::Definition() ||
        particle == geantphysics::KaonMinus::Definition() || particle == geantphysics::KaonZero::Definition() ||
        particle == geantphysics::KaonShort::Definition() || particle == geantphysics::KaonLong::Definition()) {
      // create hadronic elastic process for proton:
      //
      geantphysics::HadronicProcess *helProc = new geantphysics::ElasticScatteringProcess();
      // create the diffuse elastic model for elastic scattering
      geantphysics::HadronicFinalStateModel *diffelModel = new geantphysics::DiffuseElasticModel();
      // create the cross sections
      geantphysics::HadronicCrossSection *ggElasticXS = new geantphysics::GlauberGribovElasticXsc();

      // set min/max energies of the model
      diffelModel->SetLowEnergyUsageLimit(100.0 * geant::units::eV);
      diffelModel->SetHighEnergyUsageLimit(100.0 * geant::units::TeV);
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

} // userapplication
