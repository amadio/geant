
#include "Geant/PostStepActionStage.h"

// from geantV
#include "Geant/Propagator.h"
#include "Geant/TaskData.h"
#include "Geant/Track.h"
#include "Geant/Handler.h"

// from realphysics
#include "Geant/Material.h"
#include "Geant/MaterialCuts.h"
#include "Geant/Region.h"

#include "Geant/PhysicsProcess.h"
#include "Geant/PhysicsManagerPerParticle.h"
#include "Geant/LightTrack.h"

// handler(s)
#include "Geant/PostStepActionHandler.h"

namespace geantphysics {

PostStepActionStage::PostStepActionStage(geant::Propagator *prop)
    : geant::SimulationStage(geant::kPostStepActionStage, prop)
{
}

// base class will delete the created handlers
PostStepActionStage::~PostStepActionStage() {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PostStepActionStage::PostStepActionStage(const PostStepActionStage &other) : geant::SimulationStage(other) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PostStepActionStage &PostStepActionStage::operator=(const PostStepActionStage &other)
{
  geant::SimulationStage::operator=(other);
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
geant::SimulationStage *PostStepActionStage::Clone() const
{
  PostStepActionStage *stage = new PostStepActionStage(*this);
  return stage;
}

int PostStepActionStage::CreateHandlers()
{
  int threshold     = fPropagator->fConfig->fNperBasket;
  int physthreshold = (fPropagator->fConfig->fNvecPHY > 0) ? fPropagator->fConfig->fNvecPHY : threshold;
  // create the only one handler
  AddHandler(new PostStepActionHandler(physthreshold, fPropagator));
  // only one handler is created
  return 1;
}

// Selects tracks that have any processes, any post step processes i.e. discrete part and that limited the step
geant::Handler *PostStepActionStage::Select(geant::Track *track, geant::TaskData * /*td*/)
{
  if (track->Status() == geant::TrackStatus_t::kPhysics && track->EIndex() == 1000) {
    // these tracks should always have psorcesses active in the given region moreover should always have discrete
    // processes that limited the step (fEindex==1000)
    return fHandlers[0];
  }
  // not physics or not discrete part of limited the step
  return nullptr;
}

} // namespace geantphysics
