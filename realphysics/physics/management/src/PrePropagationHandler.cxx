
#include "PrePropagationHandler.h"

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
#include "MSCProcess.h"

#include "Geant/NavigationInterface.h"


namespace geantphysics {

PrePropagationHandler::PrePropagationHandler(int threshold, Geant::GeantPropagator *propagator)
: Geant::Handler(threshold, propagator) {}


PrePropagationHandler::~PrePropagationHandler() {}


// The PrePropagationStage will select only tracks with particles that (1) has any physics processes
// active in the given region and (2) has msc process
void PrePropagationHandler::DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td) {
  // ---
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
  assert(pManager->HasMSCProcess()); // (2)
  //
  // invoke the step limit part of the msc process
  const MSCProcess *mscProc = static_cast<const MSCProcess*>(pManager->GetMSCProcess());
  assert(mscProc);  // make sure that it is not nullptr
  // invoke the step limit method
  //  no extra geometry call anymore: geometry stage is invoked before the prepropagation stage
  //  Geant::ScalarNavInterfaceVGM::NavFindNextBoundaryMSC(*track, track->fPstep);
  mscProc->AlongStepLimitationLength(track, td);
  // --
  // copy the input track to the output
  // (secondaries should be copied as well but there is no secondary production in the msc along-step limit)
  output.AddTrack(track);
}

//______________________________________________________________________________
void PrePropagationHandler::DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td)
{
  // For the moment just loop and call scalar DoIt
  Geant::TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}


}  // namespace geantphysics
