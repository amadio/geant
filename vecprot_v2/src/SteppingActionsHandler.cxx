#include "SteppingActionsHandler.h"

#include "GeantTaskData.h"
#include "GeantVApplication.h"

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
  if (fPropagator->fStdApplication)
    fPropagator->fStdApplication->SteppingActions(*track, td);
  fPropagator->fApplication->SteppingActions(*track, td);
  
  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SteppingActionsHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector handling of stepping actions.
  
  TrackVec_t &tracks = input.Tracks();
  if (fPropagator->fStdApplication)
    fPropagator->fStdApplication->SteppingActions(tracks, td);
  fPropagator->fApplication->SteppingActions(tracks, td);

  // Copy tracks to output
#ifndef VECCORE_CUDA
  std::copy(tracks.begin(), tracks.end(), std::back_inserter(output.Tracks()));
#else
  for (auto track : tracks) output->AddTrack(track);
#endif
}

} // GEANT_IMPL_NAMESPACE
} // Geant
