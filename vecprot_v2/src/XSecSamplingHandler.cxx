#include "XSecSamplingHandler.h"

#include "PhysicsProcessOld.h"
#include "GeantTaskData.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
XSecSamplingHandler::XSecSamplingHandler(int threshold, GeantPropagator *propagator)
               : Handler(threshold, propagator)
{
// Default constructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
XSecSamplingHandler::~XSecSamplingHandler()
{
// Destructor
}  


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void XSecSamplingHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Invoke scalar xsec sampling
  // First reset step and energy deposit
  track->fStep = 0.;
  track->fEdep = 0.;
  fPropagator->fProcess->ComputeIntLen(track, td);
  // Copy to output
  output.AddTrack(track);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void XSecSamplingHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector geometry length computation. The tracks are moved into the output basket.
  TrackVec_t &tracks = input.Tracks();
  // First reset step and energy deposit
  for (auto track : tracks) {
    track->fStep = 0.;
    track->fEdep = 0.;
  }
  fPropagator->fProcess->ComputeIntLen(tracks, td);
  // Copy tracks to output
#ifndef VECCORE_CUDA
  std::copy(tracks.begin(), tracks.end(), std::back_inserter(output.Tracks()));
#else
  for (auto track : tracks) output.AddTrack(track);
#endif
}

} // GEANT_IMPL_NAMESPACE
} // Geant
