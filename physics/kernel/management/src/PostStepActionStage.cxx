
#include "PostStepActionStage.h"

// from geantV
#include "Geant/Propagator.h"
#include "Geant/TaskData.h"
#include "Geant/Track.h"
#include "Geant/Handler.h"

// from realphysics
#include "Material.h"
#include "MaterialCuts.h"
#include "Geant/Region.h"

#include "PhysicsProcess.h"
#include "PhysicsManagerPerParticle.h"
#include "LightTrack.h"

// handler(s)
#include "PostStepActionHandler.h"


namespace geantphysics {

PostStepActionStage::PostStepActionStage(geant::Propagator *prop)
: SimulationStage(geant::kPostStepActionStage, prop) {}


// base class will delete the created handlers
PostStepActionStage::~PostStepActionStage() {}


int PostStepActionStage::CreateHandlers() {
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new PostStepActionHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}


// Selects tracks that have any processes, any post step processes i.e. discrete part and that limited the step
geant::Handler* PostStepActionStage::Select(geant::Track *track, geant::TaskData * /*td*/) {
  if (track->Status() == geant::TrackStatus_t::kPhysics && track->EIndex() == 1000) {
    // these tracks should always have psorcesses active in the given region moreover should always have discrete
    // processes that limited the step (fEindex==1000)
    return fHandlers[0];
  }
  // not physics or not discrete part of limited the step
  return nullptr;
}



}  // namespace geantphysics
