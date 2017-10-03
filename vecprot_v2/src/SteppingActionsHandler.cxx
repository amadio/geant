#include "SteppingActionsHandler.h"

#include "Geant/Error.h"
#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "TrackManager.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SteppingActionsHandler::SteppingActionsHandler(int threshold, GeantPropagator *propagator)
               : Handler(threshold, propagator)
{
// Default constructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SteppingActionsHandler::~SteppingActionsHandler()
{
// Destructor
}  


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SteppingActionsHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Invoke scalar handling. Users may change the fate of the track by changing the fStage field.
  // If track made too many steps, deposit all kinetic energy and kill it
  if (track->fStatus != kNew) track->fNsteps++;
  if (track->fNsteps > fPropagator->fConfig->fNstepsKillThr) {
    Error("SteppingActions", "track %d from event %d looping -> killing it",
          track->fParticle, track->fEvent);
    track->fStatus = kKilled;
    track->Stop();
  }
  
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (fPropagator->fStdApplication)
    fPropagator->fStdApplication->SteppingActions(*track, td);
  fPropagator->fApplication->SteppingActions(*track, td);
#endif

  // The track may die at the end of the step
  if (track->fStatus == kKilled || track->fStatus == kExitingSetup) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    if (fPropagator->fStdApplication)
      fPropagator->fStdApplication->FinishTrack(*track, td);
    fPropagator->fApplication->FinishTrack(*track, td);
    fPropagator->StopTrack(track, td);
    fPropagator->fTrackMgr->ReleaseTrack(*track);
#endif
    // Dead tracks are not copied to output
    return;
  }
  
  // Update the particle location after the step
  if (track->fStatus == kBoundary)
    track->UpdateSwapPath();
  // Reset number of boundary steps
  //track->fNsteps = 0;

  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SteppingActionsHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector handling of stepping actions.
  
  TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    track->fNsteps++;
    if (track->fNsteps > fPropagator->fConfig->fNstepsKillThr) {
      Error("SteppingActions", "track %d from event %d looping -> killing it",
            track->fParticle, track->fEvent);
      track->fStatus = kKilled;
      track->Stop();
    }
  }

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (fPropagator->fStdApplication)
    fPropagator->fStdApplication->SteppingActions(tracks, td);
  fPropagator->fApplication->SteppingActions(tracks, td);
#endif

  // Copy tracks alive to output, stop the others.
  for (auto track : tracks) {
    if (track->fStatus == kKilled || track->fStatus == kExitingSetup) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
      if (fPropagator->fStdApplication)
        fPropagator->fStdApplication->FinishTrack(*track, td);
      fPropagator->fApplication->FinishTrack(*track, td);
      fPropagator->StopTrack(track, td);
      fPropagator->fTrackMgr->ReleaseTrack(*track);
#endif
    // Dead tracks are not copied to output
      continue;
    } 
       
    // Update the particle location after the step
    if (track->fStatus == kBoundary)
      track->UpdateSwapPath();
    // Reset number of boundary steps
    //track->fNsteps = 0;
    
    output.AddTrack(track);
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
