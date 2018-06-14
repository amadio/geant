
#include "Geant/PostStepActionPhysModelHandler.h"

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

PostStepActionPhysModelHandler::PostStepActionPhysModelHandler(int threshold, geant::Propagator *propagator,
                                                               int modelIdx)
    : geant::Handler(threshold, propagator)
{
  fMayBasketize = true;
  fModel        = EMModel::GetGlobalTable()[modelIdx];
}

PostStepActionPhysModelHandler::~PostStepActionPhysModelHandler()
{
}

void PostStepActionPhysModelHandler::DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td)
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
  const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
      (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
  primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
  primaryLT.SetKinE(track->Ekin());
  primaryLT.SetLogKinE(track->LogEkin());
  primaryLT.SetMass(track->Mass());
  primaryLT.SetGVcode(track->GVcode());
  primaryLT.SetDirX(track->Dx());
  primaryLT.SetDirY(track->Dy());
  primaryLT.SetDirZ(track->Dz());

  // clean the number of secondary tracks used (in PhysicsData)
  td->fPhysicsData->ClearSecondaries();

  int nSecParticles = fModel->SampleSecondaries(primaryLT, td);

  // update Track
  double newEkin = primaryLT.GetKinE();
  track->SetMass(primaryLT.GetMass());
  track->SetEkin(newEkin);
  track->SetP(std::sqrt(newEkin * (newEkin + 2.0 * track->Mass())));
  track->SetDirection(primaryLT.GetDirX(), primaryLT.GetDirY(), primaryLT.GetDirZ());
  track->SetEdep(track->Edep() + primaryLT.GetEnergyDeposit());

  if (newEkin <= 0.) {
    const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
        (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
    const Particle *particle            = Particle::GetParticleByInternalCode(track->GVcode());
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    if (primaryLT.GetTrackStatus() == LTrackStatus::kKill || pManager->GetListAtRestCandidateProcesses().size() == 0) {
      track->Kill();
      track->SetStage(geant::kSteppingActionsStage);
    } else {
      track->SetStage(geant::kAtRestActionStage);
    }
  }

  // create secondary tracks if there are any
  if (nSecParticles) {
    LightTrack *secLt = td->fPhysicsData->GetListOfSecondaries();

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
      geantTrack.SetEkin(secEkin);                                               // kinetic E of this secondary particle
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

void PostStepActionPhysModelHandler::DoItVector(geant::Track **gtracks, int N, geant::Basket &output,
                                                geant::TaskData *td)
{
  assert(N % geant::kVecLenD == 0);
  LightTrack_v &primaryLTs = td->fPhysicsData->GetPrimarySOA();

  for (int i = 0; i < N; ++i) {
    geant::Track *track = gtracks[i];

    // we will use members:
    //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
    //  fKinE         <==>  fE-fMass  // kinetic energy;  will be set to the new kinetic energy
    //  fMass         <==>  fMass     // dynamic mass of the particle
    //  fGVcode       <==>  fGVcode   // internal particle code
    //  fIntLen       <==>  fIntLen   // pre-step lambda for accounting energy loss i.e. to see if it is a delta inter.
    //  fXdir         <==>  fXdir     // direction vector x comp. will be set to the new direction x comp.
    //  fYdir         <==>  fYdir     // direction vector y comp. will be set to the new direction y comp.
    //  fZdir         <==>  fZdir     // direction vector z comp. will be set to the new direction z comp.
    const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
        (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
    primaryLTs.SetMaterialCutCoupleIndex(matCut->GetIndex(), i);
    primaryLTs.SetKinE(track->Ekin(), i);
    primaryLTs.SetLogKinE(track->LogEkin(), i);
    primaryLTs.SetMass(track->Mass(), i);
    primaryLTs.SetGVcode(track->GVcode(), i);
    primaryLTs.SetDirX(track->Dx(), i);
    primaryLTs.SetDirY(track->Dy(), i);
    primaryLTs.SetDirZ(track->Dz(), i);
    primaryLTs.SetTrackIndex(i, i);
  }
  primaryLTs.SetNtracks(N);

  // clean the number of secondary tracks used (in PhysicsData)
  td->fPhysicsData->GetSecondarySOA().ClearTracks();

  fModel->SampleSecondaries(primaryLTs, td);

  // update primary tracks
  for (int i = 0; i < N; ++i) {
    geant::Track *track = gtracks[i];
    double newEkin      = primaryLTs.GetKinE(i);
    track->SetMass(primaryLTs.GetMass(i));
    track->SetEkin(newEkin);
    track->SetP(std::sqrt(newEkin * (newEkin + 2.0 * track->Mass())));
    track->SetDirection(primaryLTs.GetDirX(i), primaryLTs.GetDirY(i), primaryLTs.GetDirZ(i));
    track->SetEdep(track->Edep() + primaryLTs.GetEnergyDeposit(i));

    //    Assert that photoelectric photons are killed
    if (track->GetPhysicsProcessIndex() == 2 && track->GVcode() == 42) {
      assert(primaryLTs.GetTrackStatus(i) == LTrackStatus::kKill);
      assert(newEkin <= 0.);
    }

    if (newEkin <= 0.) {
      const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
          (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
      const Particle *particle            = Particle::GetParticleByInternalCode(track->GVcode());
      PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
      if (primaryLTs.GetTrackStatus(i) == LTrackStatus::kKill ||
          pManager->GetListAtRestCandidateProcesses().size() == 0) {
        track->Kill();
        track->SetStage(geant::kSteppingActionsStage);
      } else {
        track->SetStage(geant::kAtRestActionStage);
      }
    }

    // copy the updated primary track to the output basket as well
    output.AddTrack(track);
  }

  //
  // create secondary tracks if there are any
  LightTrack_v &secondaryLTs = td->fPhysicsData->GetSecondarySOA();
  for (int i = 0; i < secondaryLTs.GetNtracks(); ++i) {
    auto track    = gtracks[secondaryLTs.GetTrackIndex(i)];
    int secGVcode = secondaryLTs.GetGVcode(i); // GV index of this secondary particle

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
    geantTrack.SetMass(secondaryLTs.GetMass(i));
    geantTrack.SetPosition(track->X(), track->Y(), track->Z());
    geantTrack.SetDirection(secondaryLTs.GetDirX(i), secondaryLTs.GetDirY(i), secondaryLTs.GetDirZ(i));
    double secEkin = secondaryLTs.GetKinE(i);
    geantTrack.SetP(std::sqrt(secEkin * (secEkin + 2.0 * geantTrack.Mass()))); // momentum of this secondadry particle
    geantTrack.SetEkin(secEkin);                                               // kinetic E of this secondary particle
    geantTrack.SetTime(track->Time());                                         // global time
    geantTrack.SetSafety(track->GetSafety());
    geantTrack.SetBoundary(track->Boundary());
    geantTrack.SetPath(track->Path());
    geantTrack.SetNextPath(track->Path());
    geantTrack.SetMother(track->Particle());
    geantTrack.SetPrimaryParticleIndex(track->PrimaryParticleIndex());
    // add Track
    td->fPropagator->AddTrack(geantTrack);
    output.AddTrack(&geantTrack);
  }
}

void PostStepActionPhysModelHandler::DoItScalar(geant::Track **gtracks, int N, geant::Basket &output,
                                                geant::TaskData *td)
{
  for (int i = 0; i < N; ++i) {
    geant::Track *track = gtracks[i];
    DoIt(track, output, td);
  }
}

//______________________________________________________________________________
void PostStepActionPhysModelHandler::DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td)
{
  int vectSize                                  = (input.GetNtracks() / geant::kVecLenD) * geant::kVecLenD;
  if (vectSize <= 2 * geant::kVecLenD) vectSize = 0;
  if (vectSize > 0) {
    DoItVector(input.Tracks().data(), vectSize, output, td);
  }
  int nonVectSize = input.GetNtracks() - vectSize;
  if (nonVectSize > 0) {
    DoItScalar(input.Tracks().data() + vectSize, nonVectSize, output, td);
  }
}

} // namespace geantphysics
