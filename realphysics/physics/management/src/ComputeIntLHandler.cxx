
#include "ComputeIntLHandler.h"

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

namespace geantphysics {

ComputeIntLHandler::ComputeIntLHandler(int threshold, Geant::GeantPropagator *propagator)
: Geant::Handler(threshold, propagator) {}


ComputeIntLHandler::~ComputeIntLHandler() {}

//
// The select method of the ComputeIntLStage selected only those tracks that (1) that has any process assigned i.e.
// their PhysicsManagerPerParticle object is not null in the given region and (2) they has AlongStep and/or PostStep
// processes that can limit the step.
void ComputeIntLHandler::DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td) {
  // reset step length and energy deposit
  track->SetStep(0.); // no setter for this member in GeantTrack
  track->SetEdep(0.);
  // ---
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int   particleCode         = track->GVcode();
  const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  // put some asserts here to make sure (1) that the partcile has any processes, (2) the particle has at least one
  // process with continuous or discrete parts.
  assert(pManager!=nullptr); // (1)
  assert(pManager->GetListAlongStepProcesses().size()+pManager->GetListPostStepCandidateProcesses().size()!=0);
  //
  // compute the intercation length:
  pManager->ComputeIntLen(track,td);
  // copy input track to the output
  output.AddTrack(track);
}

//______________________________________________________________________________
void ComputeIntLHandler::DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td)
{
  // For the moment just loop and call scalar DoIt
  Geant::TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}


} // namespace geantphysics
