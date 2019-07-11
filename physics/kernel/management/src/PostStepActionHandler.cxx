
#include "Geant/PostStepActionHandler.h"

// from geantV
#include "Geant/Propagator.h"
#include "Geant/TaskData.h"
#include "Geant/Track.h"
#include "Geant/Basket.h"

// from realphysics
#include "Geant/Material.h"
#include "Geant/MaterialCuts.h"
#include "Geant/Region.h"
#include "Geant/Particle.h"

#include "Geant/PhysicsProcess.h"
#include "Geant/PhysicsManagerPerParticle.h"
#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

namespace geantphysics {

PostStepActionHandler::PostStepActionHandler(int threshold, geant::Propagator *propagator)
    : geant::Handler(threshold, propagator)
{
  SetName("PostStepAction");
}

PostStepActionHandler::~PostStepActionHandler() {}

void PostStepActionHandler::DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td)
{
  // ---
  int numSecondaries = 0;
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
      (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int particleCode         = track->GVcode();
  const Particle *particle = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  // put some asserts here to make sure (1) that the partcile has any processes, (2) the particle has at least one
  // process with discrete
  assert(pManager != nullptr);                                       // (1)
  assert(pManager->GetListPostStepCandidateProcesses().size() != 0); // (2)
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
  //  pre-step mfp of the selected discrete process will also be set in pManagerPerParticle::PostStepAction!
  primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
  primaryLT.SetKinE(track->E() - track->Mass());
  primaryLT.SetLogKinE(track->LogEkin());
  primaryLT.SetMass(track->Mass());
  primaryLT.SetGVcode(track->GVcode());
  //  primaryLT.SetTrackIndex(i);
  primaryLT.SetDirX(track->Dx());
  primaryLT.SetDirY(track->Dy());
  primaryLT.SetDirZ(track->Dz());
  // obtain and set the pre-step point mfp of the selected discrete process:
  // global index of the selected physics process =  track->GetPhysicsProcessIndex();
  primaryLT.SetTotalMFP(track->GetPhysicsInteractLength(track->GetPhysicsProcessIndex()));
  //
  // clean the number of secondary tracks used (in PhysicsData)
  td->fPhysicsData->ClearSecondaries();
  //
  // invoke the PostStepAction of this particle PhysicsManagerPerParticle
  int nSecParticles = pManager->PostStepAction(primaryLT, track, td);
  //
  // update Track
  double newEkin = primaryLT.GetKinE();
  track->SetMass(primaryLT.GetMass());
  track->SetEkin(newEkin);
  track->SetP(std::sqrt(newEkin * (newEkin + 2.0 * track->Mass())));
  track->SetDirection(primaryLT.GetDirX(), primaryLT.GetDirY(), primaryLT.GetDirZ());
  track->SetEdep(track->Edep() + primaryLT.GetEnergyDeposit());
  if (newEkin <= 0.) {
    if (pManager->GetListAtRestCandidateProcesses().size() > 0 && primaryLT.GetTrackStatus() != LTrackStatus::kKill) {
      // send it to the AtRestAction stage
      track->SetStage(geant::kAtRestActionStage);
    } else {
      // kill the primary track and send the track to the last i.e. steppin-action stage
      track->Kill();
      track->SetStage(geant::kSteppingActionsStage);
    }
  }
  //
  // create secondary tracks if there are any
  if (nSecParticles) {
    // get the list of secondary tracks
    LightTrack *secLt = td->fPhysicsData->GetListOfSecondaries();
    for (int isec = 0; isec < nSecParticles; ++isec) {
      int secGVcode               = secLt[isec].GetGVcode(); // GV index of this secondary particle
      const Particle *secParticle = Particle::GetParticleByInternalCode(secGVcode);
      // get a Track geantTrack;
      geant::Track &geantTrack = td->GetNewTrack();
      // td->fNinflight++;
      // set the new track properties
      //      int t = secLt[isec].GetTrackIndex();          // parent Track index in the input Track_v
      geantTrack.SetEvent(track->Event());
      geantTrack.SetEvslot(track->EventSlot());
      geantTrack.SetGVcode(secGVcode);
      geantTrack.SetCharge(secParticle->GetPDGCharge());
      // set the index of the process (in the global process vector) that limited the step i.e. generated this secondary
      geantTrack.SetProcess(track->Process());
      geantTrack.SetStatus(geant::kNew);                 // secondary is a new track
      geantTrack.SetStage(geant::kSteppingActionsStage); // send this to the stepping action stage
      geantTrack.SetGeneration(track->GetGeneration() + 1);
      geantTrack.SetMass(secLt[isec].GetMass());
      geantTrack.SetPosition(track->X(), track->Y(), track->Z());
      geantTrack.SetDirection(secLt[isec].GetDirX(), secLt[isec].GetDirY(), secLt[isec].GetDirZ());
      double secEkin = secLt[isec].GetKinE();
      geantTrack.SetP(std::sqrt(secEkin * (secEkin + 2.0 * geantTrack.Mass()))); // momentum of this secondadry particle
      geantTrack.SetEkin(secEkin);                                               // total E of this secondary particle
      geantTrack.SetGlobalTime(track->GlobalTime());                             // global time
      geantTrack.SetSafety(track->GetSafety());
      geantTrack.SetBoundary(track->Boundary());
      geantTrack.SetPath(track->Path());
      geantTrack.SetNextPath(track->Path());
      geantTrack.SetMother(track->Particle());
      geantTrack.SetPrimaryParticleIndex(track->PrimaryParticleIndex());
      // add Track
      td->AddTrack(geantTrack);
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
void PostStepActionHandler::DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td)
{
  // For the moment just loop and call scalar DoIt
  geant::TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}

} // namespace geantphysics
