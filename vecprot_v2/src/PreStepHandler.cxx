#include "PreStepHandler.h"

#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "TrackManager.h"

namespace Geant {
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
// Invoke scalar BeginTrack user actions. All tracks arriving here should have
// the status set to Geant::kNew

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fPropagator->fApplication->BeginTrack(*track, td);
  // User may have decided to stop the track (fast simulation, ...)
#endif

  if (track->fStatus == kKilled) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    fPropagator->StopTrack(track);
    fPropagator->fTrackMgr->ReleaseTrack(*track);
#endif
    return;
  }
  
  // Set the status to "in flight"
  track->fStatus = kInFlight;
  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void PreStepHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector handling of stepping actions.
  
  TrackVec_t &tracks = input.Tracks();
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fPropagator->fApplication->BeginTrack(tracks, td);
#endif

  // Copy tracks alive to output, stop the others.
  for (auto track : tracks) {
    if (track->fStatus == kKilled) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
      fPropagator->StopTrack(track);
#endif
      continue;
    }    
    // Set the status to "in flight"
    track->fStatus = kInFlight;
    output.AddTrack(track);
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
