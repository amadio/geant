
#include "Geant/PostPropagationVectorHandler.h"

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

#include "Geant/MSCModel.h"

namespace geantphysics {

PostPropagationVectorHandler::PostPropagationVectorHandler(int threshold, geant::Propagator *propagator, int modelIdx)
    : geant::Handler(threshold, propagator)
{
  fModel = static_cast<MSCModel *>(EMModel::GetGlobalTable()[modelIdx]);
}

PostPropagationVectorHandler::~PostPropagationVectorHandler()
{
}

// The PostPropagationVectorStage will select only tracks with particles that (1) has any physics processes
// active in the given region and (2) has msc process
void PostPropagationVectorHandler::DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td)
{
  fModel->AlongStepDoIt(track, td);
  // Particles have been propagated by the geometrical step and after the msc AlongStepDoIt this geometrical step
  // is converted back to true step length. Time and number of interaction left must be updated by using this true
  // path length (that have been written into fStep by msc).
  track->IncreaseTime(track->TimeStep(track->GetStep()));
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
void PostPropagationVectorHandler::DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td)
{
  // For the moment just loop and call scalar DoIt
  geant::TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}

} // namespace geantphysics
