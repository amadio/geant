#include "ContinuousProcHandler.h"

#include "GeantTaskData.h"
#include "PhysicsProcessOld.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
ContinuousProcHandler::ContinuousProcHandler(int threshold, GeantPropagator *propagator)
               : Handler(threshold, propagator)
{
// Default constructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
ContinuousProcHandler::~ContinuousProcHandler()
{
// Destructor
}  


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ContinuousProcHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Invoke scalar along step processes
  int nextra_at_rest = 0;
#ifdef USE_REAL_PHYSICS
//      fPropagator->GetPhysicsInterface()->AlongStepAction(mat, output.GetNtracks(), output, nextra_at_rest, td);
#else
      fPropagator->Process()->Eloss(track, nextra_at_rest, output.Tracks(), td);
#endif
  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void ContinuousProcHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector geometry length computation. The tracks are moved into the output basket.

  // Copy tracks to output
  int nextra_at_rest = 0;
  TrackVec_t &tracks = input.Tracks();
#ifdef USE_REAL_PHYSICS
//      fPropagator->GetPhysicsInterface()->AlongStepAction(mat, output.GetNtracks(), output, nextra_at_rest, td);
#else
      fPropagator->Process()->Eloss(tracks, nextra_at_rest, output.Tracks(), td);
#endif
#ifndef VECCORE_CUDA
  std::copy(tracks.begin(), tracks.end(), std::back_inserter(output.Tracks()));
#else
  for (auto track : tracks) output.AddTrack(track);
#endif
}

} // GEANT_IMPL_NAMESPACE
} // Geant
