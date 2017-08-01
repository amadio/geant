
#include "PostStepActionHandler.h"

// from geantV
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "GeantTrack.h"
#include "Basket.h"

// from realphysics
#include "Material.h"
#include "MaterialCuts.h"
#include "Region.h"

#include "PhysicsProcess.h"
#include "PhysicsManagerPerParticle.h"
#include "LightTrack.h"
#include "PhysicsData.h"

namespace geantphysics {

PostStepActionHandler::PostStepActionHandler(int threshold, Geant::GeantPropagator *propagator)
: Geant::Handler(threshold, propagator) {}


PostStepActionHandler::~PostStepActionHandler() {}

void PostStepActionHandler::DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData * td) {
  // ---
  int numSecondaries = 0;
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int   particleCode         = track->GVcode();
  const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  // put some asserts here to make sure (1) that the partcile has any processes, (2) the particle has at least one
  // process with discrete
  assert(pManager!=nullptr); // (1)
  assert(pManager->GetListPostStepCandidateProcesses().size()!=0); // (2)
  //
  LightTrack primaryLT;
  // we will use members:
  //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
  //  fKinE         <==>  fE-fMass  // kinetic energy;  will be set to the new kinetic energy
  //  fMass         <==>  fMass     // dynamic mass of the particle
  //  fGVcode       <==>  fGVcode   // internal particle code
  //  fIntLen       <==>  fIntLen   // pre-step lambda for accounting energy loss i.e. to see if it is a delta inter.
  //  fXdir         <==>  fXdir     // direction vector x comp. will be set to the new direction x comp.
  //  fYdir         <==>  fYdir     // direction vector y comp. will be set to the new direction y comp.
  //  fZdir         <==>  fZdir     // direction vector z comp. will be set to the new direction z comp.
  primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
  primaryLT.SetKinE(track->E()-track->Mass());
  primaryLT.SetMass(track->Mass());
  primaryLT.SetGVcode(track->GVcode());
//  primaryLT.SetTrackIndex(i);
  primaryLT.SetDirX(track->DirX());
  primaryLT.SetDirY(track->DirY());
  primaryLT.SetDirZ(track->DirZ());
//  primaryLT.SetTotalMFP(track->GetIntLen());
  //
  // clean the number of secondary tracks used (in PhysicsData)
  td->fPhysicsData->SetNumUsedSecondaries(0);
  //
  // invoke the PostStepAction of this particle PhysicsManagerPerParticle
  int nSecParticles = pManager->PostStepAction(primaryLT, track, td);
  //
  // update GeantTrack
  double newEkin    = primaryLT.GetKinE();
  track->SetMass(primaryLT.GetMass());
  track->SetE(newEkin+track->Mass());
  track->SetP(std::sqrt(newEkin*(newEkin+2.0*track->Mass())));
  track->SetDirection(primaryLT.GetDirX(),primaryLT.GetDirY(),primaryLT.GetDirZ());
  track->SetEdep(track->Edep()+primaryLT.GetEnergyDeposit());
  if (primaryLT.GetTrackStatus()==LTrackStatus::kKill) {
    // !! later on we will send them to the AtRestActionStage !
    // kill the primary track and send the track to the last i.e. steppin-action stage
    track->Kill();
    track->SetStage(Geant::kSteppingActionsStage);
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
      Geant::GeantTrack &geantTrack = td->GetNewTrack();
      // set the new track properties
//      int t = secLt[isec].GetTrackIndex();          // parent GeantTrack index in the input GeantTrack_v
      geantTrack.SetEvent ( track->Event()     );
      geantTrack.SetEvslot( track->EventSlot()  );
      geantTrack.SetGVcode( secGVcode          );
      geantTrack.SetEindex( 0                  );
      geantTrack.SetCharge( secParticle->GetPDGCharge() );
      // set the index of the process (in the global process vector) that limited the step i.e. generated this secondary
      geantTrack.fProcess  = track->fProcess;
      geantTrack.fNsteps   = 0;
      geantTrack.fStatus   = Geant::kNew;                // secondary is a new track
      geantTrack.SetStage(Geant::kSteppingActionsStage); // send this to the stepping action stage
      geantTrack.fGeneration = track->fGeneration + 1;
      geantTrack.fMass     = secLt[isec].GetMass();
      geantTrack.SetPosition(track->PosX(),track->PosY(),track->PosZ());
      geantTrack.SetDirection(secLt[isec].GetDirX(),secLt[isec].GetDirY(),secLt[isec].GetDirZ());
      double secEkin       = secLt[isec].GetKinE();
      geantTrack.SetP(std::sqrt(secEkin*(secEkin+2.0*geantTrack.Mass()))); // momentum of this secondadry particle
      geantTrack.SetE(secEkin+geantTrack.Mass());                          // total E of this secondary particle
      geantTrack.fTime     = track->fTime; // global time
      geantTrack.fEdep     = 0.;
      geantTrack.fPstep    = 0.;
      geantTrack.fStep     = 0.;
      geantTrack.fSnext    = 0.;
      geantTrack.fSafety   = track->fSafety;
      geantTrack.fBoundary = track->fBoundary;
      geantTrack.fPending  = false;
      geantTrack.SetPath(track->Path());
      geantTrack.SetNextPath(track->Path());
      geantTrack.fMother   = track->fParticle;
      geantTrack.SetPrimaryParticleIndex(track->PrimaryParticleIndex());
      // add GeantTrack
      td->fPropagator->AddTrack(geantTrack);
      output.Tracks().push_back(&geantTrack);
      // increase the number of secondaries inserted
      ++numSecondaries;
    }
  }
  //---
  // copy the updated primary track to the output basket as well
  output.AddTrack(track);
}

//______________________________________________________________________________
void PostStepActionHandler::DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td)
{
  // For the moment just loop and call scalar DoIt
  Geant::TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}


}  // namespace geantphysics
