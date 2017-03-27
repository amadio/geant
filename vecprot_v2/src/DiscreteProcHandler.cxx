#include "DiscreteProcHandler.h"

#include "GeantTaskData.h"
#include "PhysicsProcessOld.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
DiscreteProcHandler::DiscreteProcHandler(int threshold, GeantPropagator *propagator)
               : Handler(threshold, propagator)
{
// Default constructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
DiscreteProcHandler::~DiscreteProcHandler()
{
// Destructor
}  


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void DiscreteProcHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Invoke scalar post step processes

 // Do post step actions for particles suffering a given process.
 // Surviving particles are added to the output array
#ifdef USE_REAL_PHYSICS
// fPropagator->GetPhysicsInterface()->PostStepAction(mat, nphys, output, ntotnext, td);
#else

  #ifdef USE_VECPHYS
//  fPropagator->fVectorPhysicsProcess->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);
  #endif
  // second: sample final states (based on the inf. regarding sampled
  //         target and type of interaction above), insert them into
  //         the track vector, update primary tracks;
  int ntotnext = 0;
  fPropagator->Process()->PostStepFinalStateSampling(track, ntotnext, output.Tracks(), td);
#endif

  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void DiscreteProcHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector post step handling.

  TrackVec_t &tracks = input.Tracks();
 // Do post step actions for particles suffering a given process.
 // Surviving particles are added to the output array
#ifdef USE_REAL_PHYSICS
// fPropagator->GetPhysicsInterface()->PostStepAction(mat, nphys, output, ntotnext, td);
#else

  #ifdef USE_VECPHYS
//  fPropagator->fVectorPhysicsProcess->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);
  #endif
  // second: sample final states (based on the inf. regarding sampled
  //         target and type of interaction above), insert them into
  //         the track vector, update primary tracks;
  int ntotnext = 0;
  fPropagator->Process()->PostStepFinalStateSampling(tracks, ntotnext, output.Tracks(), td);
#endif
  // Copy tracks to output
#ifndef VECCORE_CUDA
  std::copy(tracks.begin(), tracks.end(), std::back_inserter(output.Tracks()));
#else
  for (auto track : tracks) output.AddTrack(track);
#endif
}

} // GEANT_IMPL_NAMESPACE
} // Geant
