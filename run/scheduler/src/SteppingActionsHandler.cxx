#include "Geant/SteppingActionsHandler.h"

#include "Geant/Error.h"
#include "Geant/TaskData.h"
#include "Geant/UserApplication.h"
#include "Geant/TrackManager.h"

#define GV_VERBOSE_STEPPING 1

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SteppingActionsHandler::SteppingActionsHandler(int threshold, Propagator *propagator) : Handler(threshold, propagator)
{
  // Default constructor
  SetName("SteppingActions");
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SteppingActionsHandler::~SteppingActionsHandler()
{
  // Destructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SteppingActionsHandler::DoIt(Track *track, Basket &output, TaskData *td)
{
  // Invoke scalar handling. Users may change the fate of the track by changing the fStage field.
  // If track made too many steps, deposit all kinetic energy and kill it
  if (track->Status() != kNew) track->IncrementNsteps();
  if (track->GetNsteps() > fPropagator->fConfig->fNstepsKillThr) {
    Error("SteppingActions",
          "track %d from event %d looping -> killing it. Momentum = %7.4g , type = %d, parent = %d  primary= %d",
          track->Particle(), track->Event(), track->P(), track->GVcode(), track->Mother(),
          track->PrimaryParticleIndex());
    track->SetStatus(kKilled);
    track->Stop();
  }

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (fPropagator->fStdApplication) fPropagator->fStdApplication->SteppingActions(*track, td);
  fPropagator->fApplication->SteppingActions(*track, td);
#endif

#ifdef GV_VERBOSE_STEPPING
  if (fVerboseTracking) {
    // Optional 'tracking verbose' output
    track->Print("Verbose track/SteppingAction");
  }
#endif

  // The track may die at the end of the step
  if (track->Status() == kKilled || track->Status() == kExitingSetup) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    if (fPropagator->fStdApplication) fPropagator->fStdApplication->FinishTrack(*track, td);
    fPropagator->fApplication->FinishTrack(*track, td);
    fPropagator->StopTrack(track, td);
    fPropagator->fTrackMgr->ReleaseTrack(*track);
#endif
    // Dead tracks are not copied to output
    return;
  }

  // Update the particle location after the step
  if (track->Status() == kBoundary) track->UpdateSwapPath();
  // Reset number of boundary steps
  // track->fNsteps = 0;

  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SteppingActionsHandler::DoIt(Basket &input, Basket &output, TaskData *td)
{
  // Vector handling of stepping actions.

  TrackVec_t &tracks = input.Tracks();
  for (auto track : tracks) {
    track->IncrementNsteps();
    if (track->GetNsteps() > fPropagator->fConfig->fNstepsKillThr) {
      Error("SteppingActions", "track %d from event %d looping -> killing it", track->Particle(), track->Event());
      track->SetStatus(kKilled);
      track->Stop();
    }
  }

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (fPropagator->fStdApplication) fPropagator->fStdApplication->SteppingActions(tracks, td);
  fPropagator->fApplication->SteppingActions(tracks, td);
#endif

#ifdef GV_VERBOSE_STEPPING
  // Optional 'tracking verbose' output
  if (fVerboseTracking) {
    for (auto track : tracks) {
      track->Print("Verbose track/SteppingAction");
    }
  }
#endif

  // Copy tracks alive to output, stop the others.
  for (auto track : tracks) {
    if (track->Status() == kKilled || track->Status() == kExitingSetup) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
      if (fPropagator->fStdApplication) fPropagator->fStdApplication->FinishTrack(*track, td);
      fPropagator->fApplication->FinishTrack(*track, td);
      fPropagator->StopTrack(track, td);
      fPropagator->fTrackMgr->ReleaseTrack(*track);
#endif
      // Dead tracks are not copied to output
      continue;
    }

    // Update the particle location after the step
    if (track->Status() == kBoundary) track->UpdateSwapPath();
    // Reset number of boundary steps
    // track->fNsteps = 0;

    output.AddTrack(track);
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SteppingActionsHandler::SetVerbosity(bool val)
{
  fVerboseTracking = val;
#ifndef GV_VERBOSE_STEPPING
  std::cerr << "Warning: This installation has Stepping Verbosity turned OFF." << std::endl;
#endif
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
