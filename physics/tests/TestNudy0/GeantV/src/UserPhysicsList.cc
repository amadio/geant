
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
#include "Geant/TNudyEndfRecoPoint.h"
#include "Geant/NeutronNudyElasticModel.h"
#include "Geant/NeutronNudyInelasticModel.h"
#include "Geant/NeutronNudyFissionModel.h"
#include "Geant/NeutronNudyCaptureModel.h"
#include "Geant/NeutronNudyElasticXsec.h"
#include "Geant/NeutronNudyInelasticXsec.h"
#include "Geant/NeutronNudyFissionXsec.h"
#include "Geant/NeutronNudyCaptureXsec.h"

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
      if (particle == geantphysics::Neutron::Definition()) {
        // create neutron elastic process:
        //
        geantphysics::HadronicProcess *nelProc = new geantphysics::ElasticScatteringProcess();
        // create the ENDF based elastic model for elastic scattering
        geantphysics::HadronicFinalStateModel *nudyelModel = new geantphysics::NeutronNudyElasticModel();
        // create the cross sections
        geantphysics::HadronicCrossSection *nElasticXS = new geantphysics::NeutronNudyElasticXsec();
        
        // set min/max energies of the model
        nudyelModel->SetLowEnergyUsageLimit(1E-5 * geant::units::eV);
        nudyelModel->SetHighEnergyUsageLimit(20.0 * geant::units::MeV);
        // add the model to the process
        nelProc->AddModel(nudyelModel);
        // add the cross-sections to the process
        nelProc->AddCrossSection(nElasticXS);
        // add the process to the gamma particle
        AddProcessToParticle(particle, nelProc);
      }
    }
  }
} // userapplication
