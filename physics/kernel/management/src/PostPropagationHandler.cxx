
#include "Geant/PostPropagationHandler.h"

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
#include "Geant/MSCProcess.h"

namespace geantphysics {

PostPropagationHandler::PostPropagationHandler(int threshold, geant::Propagator *propagator)
    : geant::Handler(threshold, propagator)
{
  SetName("PostPropagation");
}

PostPropagationHandler::~PostPropagationHandler() {}

// The PostPropagationStage will select only tracks with particles that (1) has any physics processes
// active in the given region and (2) has msc process
void PostPropagationHandler::DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td)
{
  // ---
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts *>(
      (const_cast<vecgeom::LogicalVolume *>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int particleCode         = track->GVcode();
  const Particle *particle = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  // put some asserts here to make sure (1) that the partcile has any processes, (2) the particle has at least one
  // process with continuous parts.
  assert(pManager != nullptr);       // (1)
  assert(pManager->HasMSCProcess()); // (2)
  //
  // invoke the step limit part of the msc process
  const MSCProcess *mscProc = static_cast<const MSCProcess *>(pManager->GetMSCProcess());
  assert(mscProc); // make sure that it is not nullptr
  // invoke the along step doit method
  mscProc->AlongStepDoIt(track, td);
  // Particles have been propagated by the geometrical step and after the msc AlongStepDoIt this geometrical step
  // is converted back to true step length. Time and number of interaction left must be updated by using this true
  // path length (that have been written into fStep by msc).
  
  // POSTPONE the time step to the stepping actions (AG June 28. 2019)
  // track->IncreaseTime(track->TimeStep(track->GetStep()));
  
  // NOTE: we should save the previous step length in the Track and update this in the process or
  // in the process manager per particle only for the discrete processes BUT FOR THAT WE NEED TO SAVE the previous
  // step and we do it in the next step
  // track->fNintLen -= track->fStep/track->fIntLen;
  for (size_t i = 0; i < geant::kNumPhysicsProcess; ++i) {
    track->DecreasePhysicsNumOfInteractLengthLeft(i, track->GetStep() / track->GetPhysicsInteractLength(i));
  }
  // copy the input track to the output
  // (secondaries should be copied as well but there is no secondary production in the msc along-step-action)
  output.AddTrack(track);
}

//______________________________________________________________________________
void PostPropagationHandler::DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td)
{
  // For the moment just loop and call scalar DoIt
  geant::TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}

} // namespace geantphysics
