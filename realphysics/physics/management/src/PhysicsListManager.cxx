

#include  "PhysicsListManager.h"

// particles
#include  "Electron.h"
#include  "Positron.h"
#include  "Gamma.h"
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

// the dummy temporary region
#include  "Region.h"

#include  "PhysicsProcess.h"
#include  "PhysicsList.h"
#include  "PhysicsManagerPerParticle.h"
#include  "PhysicsParameters.h"

#include "ELossTableRegister.h" // to be able to clean
#include "ELossTableManager.h"

#include <iostream>

namespace geantphysics {

PhysicsListManager& PhysicsListManager::Instance() {
  static PhysicsListManager instance;
  return instance;
}


void PhysicsListManager::RegisterPhysicsList(PhysicsList *physlist, std::vector<bool> activeregionlist) {
  fPhysicsListVector.push_back(physlist);
  fActiveRegionMasks.push_back(activeregionlist);
}

void PhysicsListManager::RegisterPhysicsList(PhysicsList *physlist) {
  fPhysicsListVector.push_back(physlist);
}


void PhysicsListManager::CreateAllParticles() {
   Electron::Definition();
   Positron::Definition();
   Gamma::Definition();
   Proton::Definition();
   Neutron::Definition();
   PionPlus::Definition();
   PionMinus::Definition();
   PionZero::Definition();
   KaonPlus::Definition();
   KaonMinus::Definition();
   KaonZero::Definition();
   KaonShort::Definition();
   KaonLong::Definition();
   
   // get the particle table and loop over them: init ProcessManagerPerParticle vector elements to null
   // for each particle
   std::vector<Particle*> pTable = Particle::GetTheParticleTable();
   for (unsigned long i=0; i<pTable.size(); ++i) {
     Particle *particle = pTable[i];
     (particle->GetPhysicsManagerPerParticleVector()).resize(fNumOfRegions,nullptr);
   }
}


// will do the work to set the
void PhysicsListManager::BuildPhysicsLists() {

  // first create all particles and set the ProcessManagerPerParticleVector
  CreateAllParticles();

  // At this point we can do several checks later:
  //  - if there is any physics list given (error)
  //  - if only one physics list is given without specifying the list of active regions (see below)
  //  - is there any regions where non of the physics list are active  (warning: policy)
  //  - is there any regions where more than one pysics list are active (warning: policy)
  //  - did the user registered physics lists for more regions that are defined in the geometry

  // check if only one physics list was given without giving the active regions mask
  // then use the same physics list in all regions
  if (fPhysicsListVector.size()==1 && fActiveRegionMasks.size()==0) {
    std::vector<bool> *activeregionmask = new std::vector<bool>(fNumOfRegions,true);
    fActiveRegionMasks.push_back(*activeregionmask);
  } else {
    // error: more than one physics list without specifying the non-intersecting sets of active regions for them
  }


  // now we need to loop over the physics lists and call their Initialise methods to fill up the temporary physics
  // process vector in each particle, then we need to deliver that information to the appropriate element (depending on
  // the active regions) of the ProcessManagerPerParticle vector
  for (unsigned long i=0; i<fPhysicsListVector.size(); ++i) {
    // first we need to clear the temporary process vector of the particles because each physics list pushes
    // its processes per particle there
    std::vector<Particle*> pTable = Particle::GetTheParticleTable();
    for (unsigned long j=0; j<pTable.size(); ++j) {
      pTable[j]->ClearPhysicsProcessVector();
    }


    // create processes and assigne to particles according to the current physics list Initialise method
    fPhysicsListVector[i]->Initialize();

    // set the active regions mask of the PhysicsParameters obejct
    for (unsigned long j=0; j<fActiveRegionMasks[i].size(); ++j) {
      if (fActiveRegionMasks[i][j]) {
        //PhysicsParameters::GetListPhysicsParametersPerRegions()[j] = fPhysicsListVector[i]->GetPhysicsParameters();
        (fPhysicsListVector[i]->GetPhysicsParameters()->GetListActiveRegions()).push_back(true);
      } else {
        (fPhysicsListVector[i]->GetPhysicsParameters()->GetListActiveRegions()).push_back(false);
      }
    }

    for (unsigned long j=0; j<pTable.size(); ++j) {
      Particle *particle = pTable[j];
      std::vector<PhysicsProcess*> processVector = particle->GetPhysicsProcessVector();
      // if there is no any process assigned to the current particle
      if (processVector.size()==0)
        continue;

      // 1 create one PhysicsManagerPerParticle object for each partcile that the current PhysicsList has added at least
      //   one PhysicsProcess
      // 2 add all processes from the current physics list to this PhysicsManagerPerParticle
      //   and push the particle to the particle list of the phsyics process that the physics process is assigned to
      // 3 also set the active region indices both in the PhysicsProcess-es(a) and in the PhysicsManagerPerParticle(b)
      //   objects
      // 4 store the created PhysicsManagerPerParticle object pointer in a table (will be used to delete and initilize
      //   only once them)
      // 5 set pointers to regional PhysicsManagerPerParticle in the static Particle definition to point to this
      //   object where the current physics list is active
      //1
      PhysicsManagerPerParticle *pMPParticle = new PhysicsManagerPerParticle(particle,fPhysicsListVector[i]->GetPhysicsParameters());
      //2
      for (unsigned long pi=0; pi<processVector.size(); ++pi) {
        PhysicsProcess *physProc = processVector[pi];
        physProc->AddToListParticlesAssignedTo(particle);
        //3-a
        // clear this just to be sure
        (physProc->GetListActiveRegions()).clear();
        for (unsigned long k=0; k<fActiveRegionMasks[i].size(); ++k) { // loop over the regions
          (physProc->GetListActiveRegions()).push_back(fActiveRegionMasks[i][k]);
        }
        pMPParticle->AddProcess(physProc);
      }
      //3-b
      for (unsigned long k=0; k<fActiveRegionMasks[i].size(); ++k) {
        (pMPParticle->GetListActiveRegions()).push_back(fActiveRegionMasks[i][k]);
      }
      //4
      fPhysicsManagerPerParticleTable.push_back(pMPParticle);
      //5
      for (unsigned long k=0; k<fActiveRegionMasks[i].size(); ++k) {
        if (fActiveRegionMasks[i][k]) {
          (particle->GetPhysicsManagerPerParticleVector())[k] = pMPParticle;
        }
      }
    }
  }

  // loop over all Particle, all PhysicsManagerPerParticle and Initialize all
  for (unsigned long ipm=0; ipm<fPhysicsManagerPerParticleTable.size(); ++ipm) {
    fPhysicsManagerPerParticleTable[ipm]->Initialize();
  }

  // build the ELossTable-s if any
  ELossTableManager::Instance().BuildELossTables();

  // call PrepareForRun method for all ProcessManagerPerParticle:
  //  - if the particle has more than 1 kEnergyLoss process in the fAlongStepProcessVec only one will be kept
  // loop over all Particle, all PhysicsManagerPerParticle
  for (unsigned long ipm=0; ipm<fPhysicsManagerPerParticleTable.size(); ++ipm) {
    fPhysicsManagerPerParticleTable[ipm]->PrepareForRun();
  }

  //
  // will need to delete and cleear everything at the end:
  // -All created PhysicsProcess-s, Models-s, PhysicsManagerPerParticle-s all PhysicsList-s etc.
  //  ClearAll();
  // -And we should delete all Isotope-s, Element-s, Material-s, MaterialCuts-s and Particle-s
}



// at the end we can clear all physics processes and clear the physics process vector and physics manager per particle
// vector of the static Particle properties
// we also delete the physics lists added by the user and clear the local physics lists vector and active indices vector
void PhysicsListManager::ClearAll() {
  // delete all processes
  PhysicsProcess::ClearAllProcess();

  // delete all physics manager per Particle objects
  for (unsigned long i=0; i<fPhysicsManagerPerParticleTable.size(); ++i)
    delete fPhysicsManagerPerParticleTable[i];
  fPhysicsManagerPerParticleTable.clear();

  // loop over the particles and clear all vectors
  std::vector<Particle*> pTable = Particle::GetTheParticleTable();
  for (unsigned long i=0; i<pTable.size(); ++i) {
    Particle *particle = pTable[i];
    particle->ClearPhysicsProcessVector();  // objects are deleted above
    (particle->GetPhysicsManagerPerParticleVector()).clear(); // obejcts are deleted above
  }

  // delete all registred physics lists and clear the local container
  for (unsigned long i=0; i<fPhysicsListVector.size(); ++i) {
    delete fPhysicsListVector[i];
  }
  fPhysicsListVector.clear();

  // clear all active region masks local container
  for (unsigned long i=0; i<fActiveRegionMasks.size(); ++i) {
    fActiveRegionMasks[i].clear();
  }
  fActiveRegionMasks.clear();
  // delete all PhysicsParameters objects
  PhysicsParameters::Clear();
  // clear the ELossTableManager(will alos delete all ELossTable-s) and ELossTableRegister
  ELossTableManager::Instance().Clear();
  ELossTableRegister::Instance().Clear();
  // -And we should delete all Isotope-s, Element-s, Material-s, MaterialCuts-s and (Particle-s are singletones)
}



// just for testing
void PhysicsListManager::PrintAll() {
  // loop over each particle in the table;
  // get their ProcessManagerPerParticleVector where Vector is over regions
  // print all active processes at each reagion
  std::vector<Particle*> pTable = Particle::GetTheParticleTable();
  std::cout<<"\n============= Process Informations Per Region Per Particle ================ \n";
  for (int i=0; i<fNumOfRegions; ++i) {
    int regionIndx = i;
    std::cout<<"\n    ===================  Info for region with index = " << regionIndx << " ===================\n";
    std::cout<<"\n      PhysicsParameters for the region:\n";
    std::cout<<PhysicsParameters::GetPhysicsParametersForRegion(i);
    for (unsigned long j=0; j<pTable.size(); ++j) {
      Particle *particle = pTable[j];
      std::cout<< "\n      Active processes for Particle with Name = " << particle->GetName();
      PhysicsManagerPerParticle *pMPParticle = particle->GetPhysicsManagerPerParticlePerRegion(regionIndx);
      if (!pMPParticle) {
        std::cout<<"\n        There are no any active processes for this particle in this region."
                 <<"\n        So it's ProcessManagerPerParticle for this region is nullptr."
                 <<std::endl;
        continue; // go for the next particle
      }
      std::cout<<"\n        This ProcessManagerPerParticle is (also) active in regions: ";
      for(unsigned long k=0;k<pMPParticle->GetListActiveRegions().size();++k){
        if ((pMPParticle->GetListActiveRegions())[k])
          std::cout<<"    "<<k;
      }
      std::cout<<std::endl;

      std::vector<PhysicsProcess*> physProcVector = pMPParticle->GetListProcesses();
      for (unsigned long k=0; k<physProcVector.size(); ++k) {
        std::cout<<"\n        Process Name = " << physProcVector[k]->GetName();
        std::vector<bool> activeInVector = physProcVector[k]->GetListActiveRegions();
        std::cout<<"      ===> Active in regions =  ";
        for (unsigned long l=0; l<activeInVector.size(); ++l) {
           if (activeInVector[l])
             std::cout<<"   "<<l;
        }
      }
      std::cout<<std::endl;
    }
    std::cout<<std::endl;
  }
  std::cout<<"\n===========================================================================\n";
}



} // namespace geantphysics
