
#include "AlongStepActionHandler.h"

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

namespace geantphysics {

AlongStepActionHandler::AlongStepActionHandler(int threshold, Geant::GeantPropagator *propagator)
: Geant::Handler(threshold, propagator) {}


AlongStepActionHandler::~AlongStepActionHandler() {}


// The AlongStepActionStage will select only tracks with particles that (1) has any physics processes
// active in the given region and (2) has any continuous processes i.e. has along-step-action
void AlongStepActionHandler::DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td) {
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
  // process with continuous parts.
  assert(pManager!=nullptr); // (1)
  assert(pManager->GetListAlongStepProcesses().size()!=0); // (2)
  //
  // invoke the along step action for this track
  LightTrack primaryLT;
  // we will use members:
  //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
  //  fKinE         <==>  fE-fMass  // kinetic energy; will be set to the new kinetic energy
  //  fGVcode       <==>  fGVcode   // internal particle code
  //  fStepLength   <==>  fStep     // current step length
  //  fEdep         <==>  fEdep     // init to 0.0; will be set to the current energy deposit
  primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
  primaryLT.SetKinE(track->E()-track->Mass());
  primaryLT.SetMass(track->Mass());
  primaryLT.SetGVcode(track->GVcode());
  primaryLT.SetStepLength(track->GetStep());
  primaryLT.SetEnergyDeposit(0.0);
  int nSecParticles = pManager->AlongStepAction(primaryLT, td);
  // update GeantTrack
  double newEkin = primaryLT.GetKinE();
  track->SetMass(primaryLT.GetMass());
  track->SetE(newEkin+track->Mass());
  track->SetP(std::sqrt(newEkin*(newEkin+2.0*track->Mass())));
  track->SetEdep(track->Edep()+primaryLT.GetEnergyDeposit());
  if (primaryLT.GetTrackStatus()==LTrackStatus::kKill) {
    // !! later on we will send them to the AtRestActionStage !
    // kill the primary track and send the track to the last i.e. steppin-action stage
    track->Kill();
    track->SetStage(Geant::kSteppingActionsStage);
  }
  // no any secondary production so far so the next few lies does nothing
  numSecondaries += nSecParticles;
  if (nSecParticles) {
    //insert secondary tracks: nothing at the moment
  }
  // --
  // copy the input track to the output
  // (secondaries should be copied as well but there is no secondary production in the along-step-action so far)
  output.AddTrack(track);
}

//______________________________________________________________________________
void AlongStepActionHandler::DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td)
{
  // For the moment just loop and call scalar DoIt
  Geant::TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}


}  // namespace geantphysics
