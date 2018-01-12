
#include "PostPropagationStage.h"

// from geantV
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "GeantTrack.h"
#include "Handler.h"

// from realphysics
#include "Material.h"
#include "MaterialCuts.h"
#include "Region.h"

#include "PhysicsProcess.h"
#include "PhysicsManagerPerParticle.h"

#include "PostPropagationHandler.h"


namespace geantphysics {

PostPropagationStage::PostPropagationStage(Geant::GeantPropagator *prop)
: SimulationStage(Geant::kPostPropagationStage, prop) { }

// base class will delete the created handlers
PostPropagationStage::~PostPropagationStage() {}


int PostPropagationStage::CreateHandlers() {
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new PostPropagationHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}

// Selects tracks that have msc process
Geant::Handler* PostPropagationStage::Select(Geant::GeantTrack *track, Geant::GeantTaskData * /*td*/) {
  // here we will get the MaterialCuts from the LogicalVolume
  const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetMaterialCutsPtr()));
  // get the internal code of the particle
  int   particleCode         = track->GVcode();
  const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
  // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
  PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
  if (!pManager || !(pManager->HasMSCProcess())) {
    // No any physics process for the particle or no msc process
    //
    // True path length is identical to the geometrical one so update time and number of interaction left
    // here based on fSnext
    track->IncreaseTime(track->TimeStep(track->GetStep()));
    // NOTE: we should save the previous step length in the GeantTrack and update this in the process or
    // in the process manager per particle only for the discrete processes BUT FOR THAT WE NEED TO SAVE the previous
    // step and we do it in the next step
    // track->fNintLen -= track->fStep/track->fIntLen;
    for (size_t i=0; i<Geant::kNumPhysicsProcess; ++i) {
      track->DecreasePhysicsNumOfInteractLengthLeft(i, track->GetStep()/track->GetPhysicsInteractLength(i));
    }
    return nullptr;
  }
  // give back the only one handler that will perform the msc step limit phase actions
  return fHandlers[0];
}


}  // namespace geantphysics
