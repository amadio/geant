#include "LinearPropagationHandler.h"

#include "GeantTaskData.h"
#include "Geant/Error.h"
#include "Geant/NavigationInterface.h"

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
  if (track->fSnext < 1.E-8) td->fNsmall++;
  track->MakeStep(track->fSnext);
  // Update total number of steps
  td->fNsteps++;
  int nsmall = 0;

  if (track->fBoundary) {
    track->fStatus = kBoundary;
    // Find out location after boundary
    while ( IsSameLocation(*track, td) ) {
      nsmall++;
      if (nsmall > 10) {
        // Most likely a nasty overlap, some smarter action required. For now, just
        // kill the track.

        Error("LinearPropagator", "track %d from event %d stuck -> killing it",
              track->fParticle, track->fEvent);
        track->fStatus = kKilled;
        // Deposit track energy, then go directly to stepping actions
        track->Stop();
        track->SetStage(kSteppingActionsStage);
        td->fNkilled++;
        break;
      }
      track->MakeStep(1.E-3);
      td->fNpushed++;
    }
  } else {
    track->fStatus = kPhysics;
    // Update number of steps to physics
    td->fNphys++;
  }

  if (track->fSnext < 1.E-8) track->fSnext = 0;
  if (track->fSafety < 1.E-8) track->fSafety = 0;

  // Update time of flight and number of interaction lengths
//  track->fTime += track->TimeStep(track->fStep);
//  track->fNintLen -= track->fStep/track->fIntLen;

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
    if (track->fSnext < 1.E-8) td->fNsmall++;
    track->MakeStep(track->fSnext);
  }

  // Update total number of steps
  td->fNsteps += tracks.size();

  for (auto track : tracks) {
    int nsmall = 0;
    if (track->fBoundary) {
      track->fStatus = kBoundary;
      while ( IsSameLocation(*track, td) ) {
        nsmall++;
        if (nsmall > 10) {
          Error("LinearPropagator", "track %d from event %d stuck -> killing it",
                track->fParticle, track->fEvent);
          track->fStatus = kKilled;
          // Deposit track energy, then go directly to stepping actions
          track->Stop();
          track->SetStage(kSteppingActionsStage);
          td->fNkilled++;
          break;
        }
        track->MakeStep(1.E-3);
        td->fNpushed++;
      }
    } else {
      track->fStatus = kPhysics;
      // Update number of steps to physics
      td->fNphys++;
    }

    if (track->fSnext < 1.E-8) track->fSnext = 0;
    if (track->fSafety < 1.E-8) track->fSafety = 0;

    // Update time of flight and number of interaction lengths
//    track->fTime += track->TimeStep(track->fStep);
//    track->fNintLen -= track->fStep/track->fIntLen;
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
  (void)td;
  ScalarNavInterfaceTGeo::NavIsSameLocation(track, same);
#endif // USE_VECGEOM_NAVIGATOR
  if (same) return true;
  if (track.NextPath()->IsOutside())
    track.fStatus = kExitingSetup;
  return false;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
