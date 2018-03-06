
#include "Geant/PostStepActionPhysProcessHandler.h"

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
#include "Geant/EMModel.h"

namespace geantphysics {

PostStepActionPhysProcessHandler::PostStepActionPhysProcessHandler(int threshold, geant::Propagator *propagator,
                                                                   int modelIdx)
    : geant::Handler(threshold, propagator)
{
  fMayBasketize = true;
  fModel        = EMModel::GetGlobalTable()[modelIdx];
}

PostStepActionPhysProcessHandler::~PostStepActionPhysProcessHandler()
{
}

void PostStepActionPhysProcessHandler::DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td)
{
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
  primaryLT.SetMaterialCutCoupleIndex(track->GetMatCutIndex());
  primaryLT.SetKinE(track->E() - track->Mass());
  primaryLT.SetMass(track->Mass());
  primaryLT.SetGVcode(track->GVcode());
  primaryLT.SetDirX(track->Dx());
  primaryLT.SetDirY(track->Dy());
  primaryLT.SetDirZ(track->Dz());

  // clean the number of secondary tracks used (in PhysicsData)
  td->fPhysicsData->SetNumUsedSecondaries(0);

  int nSecParticles = fModel->SampleSecondaries(primaryLT, td);

  // update Track
  double newEkin = primaryLT.GetKinE();
  track->SetMass(primaryLT.GetMass());
  track->SetE(newEkin + track->Mass());
  track->SetP(std::sqrt(newEkin * (newEkin + 2.0 * track->Mass())));
  track->SetDirection(primaryLT.GetDirX(), primaryLT.GetDirY(), primaryLT.GetDirZ());
  track->SetEdep(track->Edep() + primaryLT.GetEnergyDeposit());

  if (newEkin <= 0.) {
    if (primaryLT.GetTrackStatus() == LTrackStatus::kKill || !track->HasAtRestAction()) {
      track->Kill();
      track->SetStage(geant::kSteppingActionsStage);
    } else {
      track->SetStage(geant::kAtRestActionStage);
    }
  }

  // create secondary tracks if there are any
  if (nSecParticles) {
    std::vector<LightTrack> &secLt = td->fPhysicsData->GetListOfSecondaries();

    for (int isec = 0; isec < nSecParticles; ++isec) {
      int secGVcode               = secLt[isec].GetGVcode(); // GV index of this secondary particle
      const Particle *secParticle = Particle::GetParticleByInternalCode(secGVcode);
      // get a Track geantTrack;
      geant::Track &geantTrack = td->GetNewTrack();

      geantTrack.SetEvent(track->Event());
      geantTrack.SetEvslot(track->EventSlot());
      geantTrack.SetGVcode(secGVcode);
      geantTrack.SetCharge(secParticle->GetPDGCharge());
      geantTrack.SetProcess(track->Process());
      geantTrack.SetStatus(geant::kNew);                 // secondary is a new track
      geantTrack.SetStage(geant::kSteppingActionsStage); // send this to the stepping action stage
      geantTrack.SetGeneration(track->GetGeneration() + 1);
      geantTrack.SetMass(secLt[isec].GetMass());
      geantTrack.SetPosition(track->X(), track->Y(), track->Z());
      geantTrack.SetDirection(secLt[isec].GetDirX(), secLt[isec].GetDirY(), secLt[isec].GetDirZ());
      double secEkin = secLt[isec].GetKinE();
      geantTrack.SetP(std::sqrt(secEkin * (secEkin + 2.0 * geantTrack.Mass()))); // momentum of this secondadry particle
      geantTrack.SetE(secEkin + geantTrack.Mass());                              // total E of this secondary particle
      geantTrack.SetTime(track->Time());                                         // global time
      geantTrack.SetSafety(track->GetSafety());
      geantTrack.SetBoundary(track->Boundary());
      geantTrack.SetPath(track->Path());
      geantTrack.SetNextPath(track->Path());
      geantTrack.SetMother(track->Particle());
      geantTrack.SetPrimaryParticleIndex(track->PrimaryParticleIndex());

      // add Track
      td->fPropagator->AddTrack(geantTrack);
      output.Tracks().push_back(&geantTrack);
    }
  }
  //---
  // copy the updated primary track to the output basket as well
  output.AddTrack(track);
}
//______________________________________________________________________________
void PostStepActionPhysProcessHandler::DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td)
{
  geant::TrackVec_t &gtracks = input.Tracks();

  std::vector<SecondariesFillInfo> secondFillInfo;
  secondFillInfo.reserve(gtracks.size());

  std::vector<LightTrack> primaryLTs;
  primaryLTs.reserve(gtracks.size());

  for (size_t i = 0; i < gtracks.size(); ++i) {
    geant::Track *track = gtracks[i];
    secondFillInfo.emplace_back(i, 0);

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
    primaryLT.SetMaterialCutCoupleIndex(track->GetMatCutIndex());
    primaryLT.SetKinE(track->E() - track->Mass());
    primaryLT.SetMass(track->Mass());
    primaryLT.SetGVcode(track->GVcode());
    primaryLT.SetDirX(track->Dx());
    primaryLT.SetDirY(track->Dy());
    primaryLT.SetDirZ(track->Dz());

    primaryLTs.push_back(primaryLT);
  }

  // clean the number of secondary tracks used (in PhysicsData)
  td->fPhysicsData->SetNumUsedSecondaries(0);

  fModel->SampleSecondariesVector(primaryLTs, secondFillInfo, td);

  // update primary tracks
  for (size_t i = 0; i < gtracks.size(); ++i) {
    geant::Track *track   = gtracks[i];
    LightTrack &primaryLT = primaryLTs[i];
    double newEkin        = primaryLT.GetKinE();
    track->SetMass(primaryLT.GetMass());
    track->SetE(newEkin + track->Mass());
    track->SetP(std::sqrt(newEkin * (newEkin + 2.0 * track->Mass())));
    track->SetDirection(primaryLT.GetDirX(), primaryLT.GetDirY(), primaryLT.GetDirZ());
    track->SetEdep(track->Edep() + primaryLT.GetEnergyDeposit());

    //    Assert that photoelectric photons are killed
    //    if (track->GetPhysicsProcessIndex() == 2 && track->GVcode() == 42) {
    //      assert(primaryLT.GetTrackStatus() == LTrackStatus::kKill);
    //      assert(newEkin <= 0.);
    //    }

    if (newEkin <= 0.) {
      if (primaryLT.GetTrackStatus() == LTrackStatus::kKill || !track->HasAtRestAction()) {
        track->Kill();
        track->SetStage(geant::kSteppingActionsStage);
      } else {
        track->SetStage(geant::kAtRestActionStage);
      }
    }
  }
  //
  // create secondary tracks if there are any
  int secOffset                  = 0;
  int usedSecondaries            = 0;
  std::vector<LightTrack> &secLt = td->fPhysicsData->GetListOfSecondaries();
  for (auto &secondaryInfo : secondFillInfo) {
    int nSecParticles = secondaryInfo.fNumSecondaries;
    auto track        = gtracks[secondaryInfo.fTrackId];
    // get the list of secondary tracks
    for (int isecLocal = 0; isecLocal < nSecParticles; ++isecLocal) {
      int isec                    = isecLocal + secOffset;
      int secGVcode               = secLt[isec].GetGVcode(); // GV index of this secondary particle
      const Particle *secParticle = Particle::GetParticleByInternalCode(secGVcode);
      // get a Track geantTrack;
      geant::Track &geantTrack = td->GetNewTrack();

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
      geantTrack.SetE(secEkin + geantTrack.Mass());                              // total E of this secondary particle
      geantTrack.SetTime(track->Time());                                         // global time
      geantTrack.SetSafety(track->GetSafety());
      geantTrack.SetBoundary(track->Boundary());
      geantTrack.SetPath(track->Path());
      geantTrack.SetNextPath(track->Path());
      geantTrack.SetMother(track->Particle());
      geantTrack.SetPrimaryParticleIndex(track->PrimaryParticleIndex());
      // add Track
      td->fPropagator->AddTrack(geantTrack);
      output.Tracks().push_back(&geantTrack);

      usedSecondaries++;
    }
    secOffset += nSecParticles;
  }
  assert(usedSecondaries == td->fPhysicsData->GetNumUsedSecondaries());

  for (auto track : gtracks) {
    // copy the updated primary track to the output basket as well
    output.AddTrack(track);
  }
}

} // namespace geantphysics
