#include "Geant/LinearPropagationHandler.h"

#include "Geant/TaskData.h"
#include "Geant/Error.h"
#include "Geant/NavigationInterface.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
LinearPropagationHandler::LinearPropagationHandler(int threshold, Propagator *propagator)
    : Handler(threshold, propagator)
{
  // Default constructor
  SetName("LinearPropagation");
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
LinearPropagationHandler::~LinearPropagationHandler()
{
  // Destructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void LinearPropagationHandler::DoIt(Track *track, Basket &output, TaskData *td)
{
  // Scalar geometry length computation. The track is moved into the output basket.

  // Do straight propagation to physics process or boundary
  if (track->GetSnext() < 1.E-8) td->fNsmall++;
  track->MakeStep(track->GetSnext());
  // Update total number of steps
  td->fNsteps++;
  int nsmall = 0;

  if (track->Boundary()) {
    track->SetStatus(kBoundary);
    // Find out location after boundary
    while (IsSameLocation(*track, td)) {
      nsmall++;
      if (nsmall > 10) {
        // Most likely a nasty overlap, some smarter action required. For now, just
        // kill the track.

        Error("LinearPropagator", "track %d from event %d stuck -> killing it", track->Particle(), track->Event());
        track->SetStatus(kKilled);
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
    track->SetStatus(kPhysics);
    // Update number of steps to physics
    td->fNphys++;
  }

  if (track->GetSnext() < 1.E-8) track->SetSnext(0);
  if (track->GetSafety() < 1.E-8) track->SetSafety(0);

  // Update time of flight and number of interaction lengths
  //  track->Time += track->TimeStep(track->fStep);
  //  track->fNintLen -= track->fStep/track->fIntLen;

  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void LinearPropagationHandler::DoIt(Basket &input, Basket &output, TaskData *td)
{
  // Vector geometry length computation. The tracks are moved into the output basket.
  TrackVec_t &tracks = input.Tracks();
  // This loop should autovectorize
  for (auto track : tracks) {
    if (track->GetSnext() < 1.E-8) td->fNsmall++;
    track->MakeStep(track->GetSnext());
  }

  // Update total number of steps
  td->fNsteps += tracks.size();

  for (auto track : tracks) {
    int nsmall = 0;
    if (track->Boundary()) {
      track->SetStatus(kBoundary);
      while (IsSameLocation(*track, td)) {
        nsmall++;
        if (nsmall > 10) {
          Error("LinearPropagator", "track %d from event %d stuck -> killing it", track->Particle(), track->Event());
          track->SetStatus(kKilled);
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
      track->SetStatus(kPhysics);
      // Update number of steps to physics
      td->fNphys++;
    }

    if (track->GetSnext() < 1.E-8) track->SetSnext(0);
    if (track->GetSafety() < 1.E-8) track->SetSafety(0);

    // Update time of flight and number of interaction lengths
    //    track->fTime += track->TimeStep(track->fStep);
    //    track->fNintLen -= track->fStep/track->fIntLen;
  }

// Copy tracks to output
#ifndef VECCORE_CUDA
  std::move(tracks.begin(), tracks.end(), std::back_inserter(output.Tracks()));
#else
  for (auto track : tracks)
    output.AddTrack(track);
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool LinearPropagationHandler::IsSameLocation(Track &track, TaskData *td) const
{
  // Query geometry if the location has changed for a track
  if (track.GetSafety() > 1.E-10 && track.GetSnext() > 1.E-10) {
    // Track stays in the same volume
    track.SetBoundary(false);
    return true;
  }
  bool same;
  vecgeom::NavigationState *tmpstate = td->GetPath();
  ScalarNavInterfaceVGM::NavIsSameLocation(track, same, tmpstate);
  if (same) return true;
  if (track.NextPath()->IsOutside()) track.SetStatus(kExitingSetup);
  return false;
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
