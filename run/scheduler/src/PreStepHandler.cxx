#include "PreStepHandler.h"

#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "TrackManager.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PreStepHandler::PreStepHandler(int threshold, GeantPropagator *propagator)
               : Handler(threshold, propagator)
{
// Default constructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PreStepHandler::~PreStepHandler()
{
// Destructor
}  


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void PreStepHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Invoke scalar BeginTrack user actions.

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (track->Status() == kNew)
    fPropagator->fApplication->BeginTrack(*track, td);
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
void PreStepHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
  // For the moment just loop and call scalar DoIt
  TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
