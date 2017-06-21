
#include "PhysicsProcessHandler.h"

#include "SystemOfUnits.h"

// realphysics material
#include "Material.h"
#include "Element.h"
#include "Region.h"
#include "MaterialCuts.h"

#include "PhysicsListManager.h"
#include "PhysicsProcess.h"
#include "PhysicsList.h"
#include "PhysicsList1.h"

#ifdef USE_VECPHYS
#include "VecPhysList.h"
#endif

#include "PhysicsManagerPerParticle.h"
#include "LightTrack.h"
#include "PhysicsData.h"

#include "GeantPropagator.h"
#include "GeantTaskData.h"

#include "PhysicsParameters.h"

#include <iostream>

namespace geantphysics {

PhysicsProcessHandler::PhysicsProcessHandler() {
 // nothing to do so far
}


PhysicsProcessHandler::~PhysicsProcessHandler() {
  PhysicsListManager::Instance().ClearAll();
  Material::ClearAllMaterials(); // will delete all Element and Isotope as well
  MaterialCuts::ClearAll();  // delete all MaterialCuts objects
  PhysicsData::ClearAll();   // delete all PhysicsData objects
}


void PhysicsProcessHandler::Initialize() {
  //
  // create all MaterialCuts
  MaterialCuts::CreateAll();
  //
  // print out the MaterialCuts-table
  std::cout<<MaterialCuts::GetTheMaterialCutsTable();
  //
  // print out the MaterialCuts-table
  std::cout<<Material::GetTheMaterialTable();
  //
  // set number of regions in the PhysicsListManager before we execute the user RegisterPhysicsList method(s)
  PhysicsListManager::Instance().SetNumberOfRegions(vecgeom::Region::GetNumberOfRegions());
  //
  //  Construct one user physics list and register it in the PhysicsListManager: this should be done in the application.
  //
  //  THIS IS VERY SIMILAR To Geant4 STYLE:
  //  We have only one physics list and the active region list vector is not provided: this only one physics list will
  //  be set to used(be active) in all regions automatically.
  //
      // this is what the user will need to do in their own application
#ifdef USE_VECPHYS
      PhysicsList *vecphysList = new VecPhysList("VecPhysList");
      PhysicsListManager::Instance().RegisterPhysicsList(vecphysList);
#else
      PhysicsList *physList1 = new PhysicsList1("Physics-List-1");
      PhysicsListManager::Instance().RegisterPhysicsList(physList1);
#endif
  //
  // after the user has created and registered their PhysicsList(s): build them and initialize the whole physics
  PhysicsListManager::Instance().BuildPhysicsLists();
  //
  // print out the PhysicsParameters obejct: we ha only one physics list so we have only one physics parameter object
  std::cout<<PhysicsParameters::GetThePhysicsParametersTable()[0];
  //
  // THIS is the place where we should construct as many PhysicsData objects as working threads and set their pointers
  // in the corresponding GeantTaskData objects. HOWEVER, unfortunately those GeantTaskData objects are constructedhere
  // at the moment only after the call to this Initialize() method (so we need to do it in the Application::Initialize()
  // method at the moment.
  //
  // GeantTaskData are constructed later
  //
//  for (int i=0; i<GeantPropagator::Instance()->fNthreads; ++i) {
//    GeantPropagator::Instance()->fThreadData[i]->fPhysicsData = new geantphysics::PhysicsData();
//  }
}


void PhysicsProcessHandler::ComputeIntLen(Material_t * /*mat*/, int ntracks, GeantTrack_v &tracks, double * /*lengths*/,
                                          GeantTaskData *td) {
  for (int i=0; i<ntracks; ++i) {
    const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(tracks.GetVolume(i))->GetMaterialCutsPtr()));
    int   particleCode         = tracks.fGVcodeV[i];
    const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
    // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    if (!pManager) {
      tracks.fPstepV[i]  = PhysicsProcess::GetAVeryLargeValue();
      tracks.fIntLenV[i] = 1.0;
      continue; // no physics limit
    }
    // check if the partcile has anything along step and if yes get the minimum of their along-step limit
    // create a LightTrack
    if (pManager->GetListAlongStepProcesses().size()+pManager->GetListPostStepCandidateProcesses().size()>0) {
       LightTrack primaryLT;
       //LightTrack *lt = &(td->fPhysicsData->GetListOfSecondaries()[0]);
       // we will use members:
       //  fNintLen      <==>  fNintLen  // number of interaction left
       //  fTargetZ      <==>  fEindex   // will be set to flag if disc. or cont. step limit won
       //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
       //  fKinE         <==>  fE-fMass  // kinetic energy
       //  fGVcode       <==>  fGVcode   // internal particle code
       //  fIntLen       <==>  fIntLen   // will be set to store the current inverse total lambda
       //  fStepLength   <==>  fPstep    // will be set to store the physics step limit
       primaryLT.SetNumOfInteractionLegthLeft(tracks.fNintLenV[i]);
       primaryLT.SetKinE(tracks.fEV[i]-tracks.fMassV[i]);
       primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
       primaryLT.SetGVcode(tracks.fGVcodeV[i]);
       pManager->ComputeIntLen(primaryLT,td);
       // update GeantTrack
       tracks.fPstepV[i]   = primaryLT.GetStepLength();
       tracks.fIntLenV[i]  = primaryLT.GetTotalMFP();
       tracks.fEindexV[i]  = primaryLT.GetTargetZ();
       if (tracks.fNintLenV[i]<=0.0) { // was resampled
         tracks.fNintLenV[i] = primaryLT.GetNumOfInteractionLegthLeft();
       }
    }
  }
}


void PhysicsProcessHandler::AlongStepAction(Material_t * /*mat*/, int ntracks, GeantTrack_v &tracks, int &nout,
                                            GeantTaskData *td) {
  int numSecondaries = 0;
  for (int i=0; i<ntracks; ++i) {
    const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(tracks.GetVolume(i))->GetMaterialCutsPtr()));
    int particleCode           = tracks.fGVcodeV[i];
    const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
    // check if the partcile has anything along step
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    if (!pManager) {
      continue; // no physics of this particle in this region
    }
    if (pManager->GetListAlongStepProcesses().size()>0) {
       LightTrack primaryLT;
       // we will use members:
       //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
       //  fKinE         <==>  fE-fMass  // kinetic energy; will be set to the new kinetic energy
       //  fGVcode       <==>  fGVcode   // internal particle code
       //  fStepLength   <==>  fStep     // current step length
       //  fEdep         <==>  fEdep     // init to 0.0; will be set to the current energy deposit
       primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
       primaryLT.SetKinE(tracks.fEV[i]-tracks.fMassV[i]);
       primaryLT.SetGVcode(tracks.fGVcodeV[i]);
       primaryLT.SetStepLength(tracks.fStepV[i]);
       primaryLT.SetEnergyDeposit(0.0);
       //primaryLT.SetTrackIndex(i);
       int nSecParticles = pManager->AlongStepAction(primaryLT, td);
       // update GeantTrack
       double newEkin = primaryLT.GetKinE();
       tracks.fEV[i]  = newEkin+tracks.fMassV[i];
       tracks.fPV[i]  = std::sqrt(newEkin*(newEkin+2.0*tracks.fMassV[i]));
       tracks.fEdepV[i] += primaryLT.GetEnergyDeposit();
       if (primaryLT.GetTrackStatus()==LTrackStatus::kKill) {
         tracks.fStatusV[i] = Geant::kKilled;
       }
       numSecondaries += nSecParticles;
       if (nSecParticles) {
         //insert secondary tracks
       }
    }
  }
  nout = numSecondaries; // number of inserted GeantTracks
}


void PhysicsProcessHandler::PostStepAction(Material_t * /*mat*/, int ntracks, GeantTrack_v &tracks, int &nout,
                                           GeantTaskData *td) {
  int numSecondaries = 0;
  for (int i=0; i<ntracks; ++i) {
    const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(tracks.GetVolume(i))->GetMaterialCutsPtr()));
    int particleCode           = tracks.fGVcodeV[i];
    const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
    // check if the partcile has anything along step
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    if (!pManager) {
      continue; // no physics of this particle in this reagion
    }
    if (pManager->GetListPostStepCandidateProcesses().size()>0) {
      LightTrack primaryLT;
      // we will use members:
      //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
      //  fKinE         <==>  fE-fMass  // kinetic energy;  will be set to the new kinetic energy
      //  fGVcode       <==>  fGVcode   // internal particle code
      //  fIntLen       <==>  fIntLen   // pre-step lambda for accounting energy loss i.e. to see if it is a delta inter.
      //  fXdir         <==>  fXdir     // direction vector x comp. will be set to the new direction x comp.
      //  fYdir         <==>  fYdir     // direction vector y comp. will be set to the new direction y comp.
      //  fZdir         <==>  fZdir     // direction vector z comp. will be set to the new direction z comp.
      primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
      primaryLT.SetKinE(tracks.fEV[i]-tracks.fMassV[i]);
      primaryLT.SetGVcode(tracks.fGVcodeV[i]);
      primaryLT.SetTrackIndex(i);
      primaryLT.SetDirX(tracks.fXdirV[i]);
      primaryLT.SetDirY(tracks.fYdirV[i]);
      primaryLT.SetDirZ(tracks.fZdirV[i]);
      primaryLT.SetTotalMFP(tracks.fIntLenV[i]);
      //
      // clean the number of secondary tracks used (in PhysicsData)
      td->fPhysicsData->SetNumUsedSecondaries(0);
      //
      // invoke the PostStepAction of this particle PhysicsManagerPerParticle
      int nSecParticles = pManager->PostStepAction(primaryLT, td);
      //
      // update GeantTrack
      double newEkin    = primaryLT.GetKinE();
      tracks.fEV[i]     = newEkin+tracks.fMassV[i];
      tracks.fPV[i]     = std::sqrt(newEkin*(newEkin+2.0*tracks.fMassV[i]));
      tracks.fXdirV[i]  = primaryLT.GetDirX();
      tracks.fYdirV[i]  = primaryLT.GetDirY();
      tracks.fZdirV[i]  = primaryLT.GetDirZ();
      tracks.fEdepV[i] += primaryLT.GetEnergyDeposit();
      if (primaryLT.GetTrackStatus()==LTrackStatus::kKill) {
        tracks.fStatusV[i] = Geant::kKilled;
      }
      //
      // create secondary tracks if there are any
      if (nSecParticles) {
        // get the list of secondary tracks
        std::vector<LightTrack> &secLt = td->fPhysicsData->GetListOfSecondaries();
        for (int isec=0; isec<nSecParticles; ++isec) {
          int   secGVcode = secLt[isec].GetGVcode(); // GV index of this secondary particle
          const Particle *secParticle = Particle::GetParticleByInternalCode(secGVcode);
          // get a GeantTrack geantTrack;
          Geant::GeantTrack &geantTrack = td->GetTrack();
          // set the new track properties
          int t = secLt[isec].GetTrackIndex();          // parent GeantTrack index in the input GeantTrack_v
          geantTrack.fEvent    = tracks.fEventV[t];
          geantTrack.fEvslot   = tracks.fEvslotV[t];
          geantTrack.fGVcode   = secGVcode;
          geantTrack.fEindex   = 0;
          geantTrack.fCharge   = secParticle->GetPDGCharge();
          geantTrack.fProcess  = 0;
          geantTrack.fNsteps   = 0;
          geantTrack.fStatus   = Geant::kNew;           // this secondary is a new track
          geantTrack.fMass     = secParticle->GetPDGMass();
          geantTrack.fXpos     = tracks.fXposV[t]; // rx of this particle (same as parent)
          geantTrack.fYpos     = tracks.fYposV[t]; // ry of this particle (same as parent)
          geantTrack.fZpos     = tracks.fZposV[t]; // rz of this particle (same as parent)
          geantTrack.fXdir     = secLt[isec].GetDirX();     // dirx of this particle (before transform.)
          geantTrack.fYdir     = secLt[isec].GetDirY();     // diry of this particle before transform.)
          geantTrack.fZdir     = secLt[isec].GetDirZ();     // dirz of this particle before transform.)
          double secEkin       = secLt[isec].GetKinE();
          geantTrack.fP        = std::sqrt(secEkin*(secEkin-2.0*geantTrack.fMass)); // momentum of this secondadry particle
          geantTrack.fE        = secEkin+geantTrack.fMass;                          // total E of this secondary particle
          geantTrack.fTime     = tracks.fTimeV[t]; // global time
          geantTrack.fEdep     = 0.;
          geantTrack.fPstep    = 0.;
          geantTrack.fStep     = 0.;
          geantTrack.fSnext    = 0.;
          geantTrack.fSafety   = tracks.fSafetyV[t];
          geantTrack.fBoundary = tracks.fBoundaryV[t];
          geantTrack.fPending  = false;
          *geantTrack.fPath    = *tracks.fPathV[t];
          *geantTrack.fNextpath= *tracks.fPathV[t];
          geantTrack.fMother   = tracks.fParticleV[t];
          // add GeantTrack
          td->fPropagator->AddTrack(geantTrack);
          tracks.AddTrack(geantTrack);
          // increase the number of secondaries inserted
          ++numSecondaries;
        }
      }
    }
  }
  // set the number of secondaries insterted
  nout = numSecondaries;
}


void PhysicsProcessHandler::AtRestAction(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks*/,
                                         int & /*nout*/, GeantTaskData * /*td*/) {
  //
  // IMPLEMENTATION NEEDED!
  //
}

} // namespace geantphysics
