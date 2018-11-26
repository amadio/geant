#include "Geant/PreStepHandler.h"

#include "Geant/TaskData.h"
#include "Geant/UserApplication.h"
#include "Geant/TrackManager.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PreStepHandler::PreStepHandler(int threshold, Propagator *propagator) : Handler(threshold, propagator)
{
  // Default constructor
  SetName("PreStep");
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PreStepHandler::~PreStepHandler()
{
  // Destructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void PreStepHandler::DoIt(Track *track, Basket &output, TaskData *td)
{
  // Invoke scalar BeginTrack user actions.

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (track->Status() == kNew) fPropagator->fApplication->BeginTrack(*track, td);
// User may have decided to stop the track (fast simulation, ...)
#endif

  if (track->Status() == kKilled) {
    // The track has to be actually killed by the stepping actions
    track->SetStage(kSteppingActionsStage);
  } else {
    // Set the status to "in flight"
    track->SetStatus(kInFlight);
    track->SetPrePropagationDone(false);
  }
  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void PreStepHandler::DoIt(Basket &input, Basket &output, TaskData *td)
{
  // For the moment just loop and call scalar DoIt
  TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
