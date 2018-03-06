
#include "Geant/PostStepActionDeltaIntHandler.h"

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

PostStepActionDeltaIntHandler::PostStepActionDeltaIntHandler(int threshold, geant::Propagator *propagator)
    : geant::Handler(threshold, propagator)
{
  fMayBasketize = true;
}

PostStepActionDeltaIntHandler::~PostStepActionDeltaIntHandler()
{
}

void PostStepActionDeltaIntHandler::DoIt(geant::Track *track, geant::Basket &output, geant::TaskData * /* td */)
{
  output.AddTrack(track);
}
//______________________________________________________________________________
void PostStepActionDeltaIntHandler::DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData * /* td */)
{
  geant::TrackVec_t &gtracks = input.Tracks();
  for (auto track : gtracks) {
    output.AddTrack(track);
  }
}

} // namespace geantphysics
