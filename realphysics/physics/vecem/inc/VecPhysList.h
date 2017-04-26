#ifndef VECPHYS_LIST
#define VECPHYS_LIST

#include "SystemOfUnits.h"

#include "PhysicsList.h"
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

#include "GUGammaComptonProcess.h"
#include "GUKleinNishinaComptonModel.h"
#include "GUGammaPhotoElectricProcess.h"

#include "GUSauterGavrilaModel.h"
#include "GUGammaConversionProcess.h"
#include "GUBetheHeitlerConversionModel.h"

namespace geantphysics {

class VecPhysList : public PhysicsList {

public:
  VecPhysList(const std::string &name) : PhysicsList(name) {}

  virtual void Initialize() {
    // get the partcile table and loop over it
    std::vector<Particle*> pTable = Particle::GetTheParticleTable();
    for (unsigned int i=0; i<pTable.size(); ++i) {
      Particle *particle = pTable[i];

      if (particle==Gamma::Definition()) {
        //1. create compton process with the Klein-Nishina model:

        EMPhysicsProcess *gCompProc = new GUGammaComptonProcess("gCompton");
        EMModel          *gKNModel  = new GUKleinNishinaComptonModel(); // true);

        // set min/max energies of the model
        gKNModel->SetLowEnergyUsageLimit (1.0*geant::keV);
        gKNModel->SetHighEnergyUsageLimit(1.0*geant::TeV);

        // add the model to the process
        gCompProc->AddModel(gKNModel);

        // add the process to the gamma particle
        std::cout<< " Adding Compton to gamma."<<std::endl;
        AddProcessToPartcile(particle, gCompProc);
        std::cout<< " Adding Compton to gamma - done."<<std::endl;
        
        //2. create photo-electri process with the SauterGavrila angular distribution:
        EMPhysicsProcess *gPhotoElecProc = new GUGammaPhotoElectricProcess("gPhotoElectic");
        EMModel          *gSGModel  = new GUSauterGavrilaModel(true);

        // set min/max energies of the model
        gSGModel->SetLowEnergyUsageLimit (1.0*geant::keV);
        gSGModel->SetHighEnergyUsageLimit(1.0*geant::GeV);

        // add the model to the process
        gPhotoElecProc->AddModel(gSGModel);

        // add the process to the gamma particle
        std::cout<< " Adding Photo-Electric to gamma." << std::endl;
        AddProcessToPartcile(particle, gPhotoElecProc);
        std::cout<< " Adding Photo-Electric to gamma - done." << std::endl;

        //3. create the conversion process with the Bethe-Heitler model
        EMPhysicsProcess *gConvProc = new GUGammaConversionProcess("gConversion");
        EMModel          *gBHModel  = new GUBetheHeitlerConversionModel();

        // set min/max energies of the model
        gBHModel->SetLowEnergyUsageLimit (2.0*geant::kElectronMassC2);
        gBHModel->SetHighEnergyUsageLimit(1.0*geant::TeV);

        // add the model to the process
        gConvProc->AddModel(gBHModel);

        // add the process to the gamma particle
        std::cout<< " Adding Conversion to gamma."<<std::endl;
        AddProcessToPartcile(particle,gConvProc );
        std::cout<< " Adding Conversion to gamma - done."<<std::endl;                
      }

      if (particle==Electron::Definition()) {
        //std::cout<<"  ELECTRON" <<std::endl;
        //
        // create ionization process for e- with 1 model:
       //
        EMPhysicsProcess *eIoniProc = new ElectronIonizationProcess("eIoniFore-");
        // create the Moller-Bhabha model for ionization i.e. for e- + e- -> e- + e- intercation
        EMModel          *eMBModel  = new MollerBhabhaIonizationModel(true);
        // set min/max energies of the model
        eMBModel->SetLowEnergyUsageLimit ( 1.0*geant::keV);
        eMBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
        // add the model to the process
        eIoniProc->AddModel(eMBModel);
        //
        // add the process to the e- particle
        AddProcessToPartcile(particle, eIoniProc);
        //
        // create bremsstrahlung process for e- with 2 models:
        //
        EMPhysicsProcess *eBremProc = new ElectronBremsstrahlungProcess("eBremFore-");
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
        eRelBModel->SetLowEnergyUsageLimit ( 1.0*geant::GeV);
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
        EMPhysicsProcess *eIoniProc = new ElectronIonizationProcess("eIoniFore+");
        // create the Moller-Bhabha model for ionization i.e. for e+ + e- -> e+ + e- intercation
        EMModel          *eMBModel  = new MollerBhabhaIonizationModel(false);
        // set min/max energies of the model
        eMBModel->SetLowEnergyUsageLimit ( 1.0*geant::keV);
        eMBModel->SetHighEnergyUsageLimit(10.0*geant::TeV);
        // add the model to the process
        eIoniProc->AddModel(eMBModel);
        // add the process to the e+ particle
        AddProcessToPartcile(particle, eIoniProc);
        //
        // create bremsstrahlung process for e+ with 2 models:
        //
        EMPhysicsProcess *eBremProc = new ElectronBremsstrahlungProcess("eBremFore+");
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
        eRelBModel->SetLowEnergyUsageLimit ( 1.0*geant::GeV);
        eRelBModel->SetHighEnergyUsageLimit(10.0*geant::TeV);
        // add this model to the process
        eBremProc->AddModel(eRelBModel);
        //
        // add the process to the e+ particle
        AddProcessToPartcile(particle, eBremProc);
      }
    }
  }
};

} // namespace geantphysics
#endif // VECPHYS_LIST
