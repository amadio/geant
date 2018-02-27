
#include "Geant/AtRestActionHandler.h"

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

AtRestActionHandler::AtRestActionHandler(int threshold, geant::Propagator *propagator)
: geant::Handler(threshold, propagator) {}


AtRestActionHandler::~AtRestActionHandler() {}

void AtRestActionHandler::DoIt(geant::Track *track, geant::Basket& output, geant::TaskData * td) {
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
  // process with at-rest part
  // NOTE: kinetic energy should be <= 0 that was checked at the AtRestActionStage
  assert(pManager!=nullptr); // (1)
  assert(pManager->GetListAtRestCandidateProcesses().size()!=0); // (2)
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
///  the primary track is stopped i.e. zero kinetic energy
  primaryLT.SetKinE(0.);
  primaryLT.SetMass(track->Mass());
  primaryLT.SetGVcode(track->GVcode());
//  primaryLT.SetTrackIndex(i);
///  the primary track is stopped i.e. zero kinetic energy
///  primaryLT.SetDirX(track->Dx());
///  primaryLT.SetDirY(track->Dy());
///  primaryLT.SetDirZ(track->Dz());
//  primaryLT.SetTotalMFP(track->GetIntLen());
  //
  // clean the number of secondary tracks used (in PhysicsData)
  td->fPhysicsData->SetNumUsedSecondaries(0);
  //
  // invoke the AtRestAction of this particle PhysicsManagerPerParticle
  int nSecParticles = pManager->AtRestAction(primaryLT, track, td);
  //
  // update Track: the primary track was stopped before the interaction (i.e. zero kinetic energy) and killed
  //                    in the interaction
///  double newEkin    = primaryLT.GetKinE();
///  track->SetMass(primaryLT.GetMass());
///  track->SetE(newEkin+track->Mass());
  track->SetE(track->Mass());
///  track->SetP(std::sqrt(newEkin*(newEkin+2.0*track->Mass())));
  track->SetP(0.);
///  track->SetDirection(primaryLT.GetDirX(),primaryLT.GetDirY(),primaryLT.GetDirZ());
  track->SetEdep(track->Edep()+primaryLT.GetEnergyDeposit());
//  if (primaryLT.GetTrackStatus()==LTrackStatus::kKill) { // as it should always be !!!
    track->Kill();
    track->SetStage(geant::kSteppingActionsStage);
//  }
  //
  // create secondary tracks if there are any
  if (nSecParticles) {
    // get the list of secondary tracks
    std::vector<LightTrack> &secLt = td->fPhysicsData->GetListOfSecondaries();
    for (int isec=0; isec<nSecParticles; ++isec) {
      int   secGVcode = secLt[isec].GetGVcode(); // GV index of this secondary particle
      const Particle *secParticle = Particle::GetParticleByInternalCode(secGVcode);
      // get a Track geantTrack;
      geant::Track &geantTrack = td->GetNewTrack();
      // set the new track properties
//      int t = secLt[isec].GetTrackIndex();          // parent Track index in the input Track_v
      geantTrack.SetEvent (track->Event());
      geantTrack.SetEvslot(track->EventSlot());
      geantTrack.SetGVcode(secGVcode);
      geantTrack.SetCharge(secParticle->GetPDGCharge());
      // set the index of the process (in the global process vector) that limited the step i.e. generated this secondary
      geantTrack.SetProcess(track->Process());
      geantTrack.SetStatus(geant::kNew);                // secondary is a new track
      geantTrack.SetStage(geant::kSteppingActionsStage); // send this to the stepping action stage
      geantTrack.SetGeneration( track->GetGeneration() + 1);
      geantTrack.SetMass(secLt[isec].GetMass());
      geantTrack.SetPosition(track->X(),track->Y(),track->Z());
      geantTrack.SetDirection(secLt[isec].GetDirX(),secLt[isec].GetDirY(),secLt[isec].GetDirZ());
      double secEkin       = secLt[isec].GetKinE();
      geantTrack.SetP(std::sqrt(secEkin*(secEkin+2.0*geantTrack.Mass()))); // momentum of this secondadry particle
      geantTrack.SetE(secEkin+geantTrack.Mass());                          // total E of this secondary particle
      geantTrack.SetTime(track->Time() ); // global time
      geantTrack.SetSafety(track->GetSafety());
      geantTrack.SetBoundary(track->Boundary());
      geantTrack.SetPath(track->Path());
      geantTrack.SetNextPath(track->Path());
      geantTrack.SetMother(track->Particle());
      geantTrack.SetPrimaryParticleIndex(track->PrimaryParticleIndex());
      // add Track
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
void AtRestActionHandler::DoIt(geant::Basket &input, geant::Basket& output, geant::TaskData *td)
{
  // For the moment just loop and call scalar DoIt
  geant::TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}


}  // namespace geantphysics
