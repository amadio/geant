#include "LinearPropagationHandler.h"

#include "GeantTaskData.h"
#include "Geant/Error.h"
#include "ScalarNavInterfaceVGM.h"
#include "ScalarNavInterfaceTGeo.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
LinearPropagationHandler::LinearPropagationHandler(int threshold, GeantPropagator *propagator)
               : Handler(threshold, propagator)
{
// Default constructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
LinearPropagationHandler::~LinearPropagationHandler()
{
// Destructor
}  


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void LinearPropagationHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Scalar geometry length computation. The track is moved into the output basket.

  // Do straight propagation to physics process or boundary
  track->fPstep -= track->fSnext;
  track->fStep += track->fSnext;
  track->fSafety -= track->fSnext;
  if (track->fSafety < 0.)
    track->fSafety = 0;
  track->fXpos += track->fSnext * track->fXdir;
  track->fYpos += track->fSnext * track->fYdir;
  track->fZpos += track->fSnext * track->fZdir;
  // Update total number of steps
  td->fNsteps++;
  track->SetStage(kContinuousProcStage);

  if (track->fBoundary) {
    track->fStatus = kBoundary;
    track->fNsteps++;
    // Find out location after boundary
    while ( IsSameLocation(*track, td) ) {
      track->fSnext = 1.E-7;
      DoIt(track, output, td);
    }
  } else {
    track->fStatus = kPhysics;
    // Update number of steps to physics
    td->fNphys++;
  }

  if (track->fSnext < 1.E-8) td->fNsmall++;
  track->fSnext = 0;
  
  // Update time of flight and number of interaction lengths
  track->fTime += track->TimeStep(track->fStep);
  track->fNintLen -= track->fStep/track->fIntLen;

  // Set follow-up stage to be ContinuousProcStage

  // Kill loopers
  if (track->fNsteps > fPropagator->fConfig->fNstepsKillThr) {
    Error("TransportTracks", "track %d from event %d seems to be stuck -> killing it",
          track->fParticle, track->fEvent);
    track->fStatus = kKilled;
    track->SetStage(kSteppingActionsStage);
  }

  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void LinearPropagationHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector geometry length computation. The tracks are moved into the output basket.
  TrackVec_t &tracks = input.Tracks();  
  // This loop should autovectorize
  for (auto track : tracks) {
    track->fPstep -= track->fSnext;
    track->fStep += track->fSnext;
    track->fSafety -= track->fSnext;
    if (track->fSafety < 0.)
      track->fSafety = 0;
    track->fXpos += track->fSnext * track->fXdir;
    track->fYpos += track->fSnext * track->fYdir;
    track->fZpos += track->fSnext * track->fZdir;
  }

  // Update total number of steps
  td->fNsteps += tracks.size();

  for (auto track : tracks) {
    // Set follow-up stage to be ContinuousProcStage
    track->SetStage(kContinuousProcStage);
    if (track->fBoundary) {
      track->fStatus = kBoundary;
      track->fNsteps++;
      while ( IsSameLocation(*track, td) ) {
        track->fSnext = 1.E-7;
        DoIt(track, output, td);
        return;
      }
    } else {
      track->fStatus = kPhysics;
      // Update number of steps to physics
      td->fNphys++;
    }
    if (track->fSnext < 1.E-8) td->fNsmall++;
    track->fSnext = 0;    
    // Update time of flight and number of interaction lengths
    track->fTime += track->TimeStep(track->fStep);
    track->fNintLen -= track->fStep/track->fIntLen;  
    // Kill loopers
    if (track->fNsteps > fPropagator->fConfig->fNstepsKillThr) {
      Error("TransportTracks", "track %d from event %d seems to be stuck -> killing it",
            track->fParticle, track->fEvent);
      track->fStatus = kKilled;
      track->SetStage(kSteppingActionsStage);
    }
  }
  

  // Copy tracks to output
#ifndef VECCORE_CUDA
  std::move(tracks.begin(), tracks.end(), std::back_inserter(output.Tracks()));
#else
  for (auto track : tracks) output.AddTrack(track);
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool LinearPropagationHandler::IsSameLocation(GeantTrack &track, GeantTaskData *td) {
// Query geometry if the location has changed for a track
  if (track.fSafety > 1.E-10 && track.fSnext > 1.E-10) {
    // Track stays in the same volume
    track.fBoundary = false;
    return true;
  }  
  bool same;
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::NavigationState *tmpstate = td->GetPath();
  ScalarNavInterfaceVGM::NavIsSameLocation(track, same, tmpstate);
#else
// ROOT navigation
  ScalarNavInterfaceTGeo::NavIsSameLocation(track, same);
#endif // USE_VECGEOM_NAVIGATOR
  if (same) return true;
  if (track.fNextpath->IsOutside())
    track.fStatus = kExitingSetup;
  return false;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
