#include "GeomQueryHandler.h"
#include "Basket.h"
#include "GeantTaskData.h"
#include "GeantTrackGeo.h"
#include "GeantRunManager.h"
#include "VBconnector.h"

#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
#include "GeantNuma.h"
#endif

#ifdef USE_VECGEOM_NAVIGATOR
#include "ScalarNavInterfaceVGM.h"
#include "VectorNavInterface.h"
#else
#include "ScalarNavInterfaceTGeo.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeomQueryHandler::GeomQueryHandler(Volume_t *vol, int threshold, GeantPropagator *propagator, int index)
               : Handler(threshold, propagator), fVolume(vol), fIndex(index)
{
// Default constructor
  assert(vol && "GeomQueryHandler: A valid volume pointer has to be provided");
  ConnectToVolume();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeomQueryHandler::~GeomQueryHandler()
{
// Destructor
  DisconnectVolume();
}  

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeomQueryHandler::ConnectToVolume()
{
  VBconnector *connector = new VBconnector(fIndex);
#ifdef USE_VECGEOM_NAVIGATOR
    fVolume->SetBasketManagerPtr(connector);
#else
    fVolume->SetFWExtension(connector);
#endif
  
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeomQueryHandler::ActivateBasketizing(bool flag)
{
// Special basketizing in case of logical volumes
  if (fActive == flag) return;
  fActive = flag;
  int basket_size = fThreshold;
  // Set a 'compromise' size for the basketizer buffer for all geometry handlers
  int buffer_size = 1 << 10; // This makes ~16 MBytes of basketizer buffers for 4K volumes
  assert(fThreshold < 512 && fThreshold < buffer_size);
  // Create basketizer the first time the handler is activated
  if (fActive && !fBasketizer) {
#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
    if (GetNode() < 0) {
      fBasketizer = new basketizer_t(buffer_size, basket_size);
    } else {
      int basketizer_size = basketizer_t::SizeofInstance(buffer_size);
      fBasketizer = basketizer_t::MakeInstanceAt(
        NumaUtils::NumaAlignedMalloc(basketizer_size, GetNode(), 64), buffer_size, basket_size);
    }
#else
    fBasketizer = new basketizer_t(buffer_size, basket_size);
#endif
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeomQueryHandler::DisconnectVolume()
{
  if (!fVolume) return;
  VBconnector *connector = nullptr;
#ifdef USE_VECGEOM_NAVIGATOR
    connector = static_cast<VBconnector*>(fVolume->GetBasketManagerPtr());
    fVolume->SetBasketManagerPtr(nullptr);
#else
    connector = static_cast<VBconnector*>(fVolume->GetFWExtension());
    fVolume->SetFWExtension(nullptr);
#endif
    delete connector;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeomQueryHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Scalar geometry length computation. The track is moved into the output basket.

  // Below we should call just finding the next boundary. Relocation should
  // be handled separately
#ifdef USE_VECGEOM_NAVIGATOR
  ScalarNavInterfaceVGM::NavFindNextBoundary(*track);
#else
// ROOT geometry
  ScalarNavInterfaceTGeo::NavFindNextBoundary(*track);
#endif // USE_VECGEOM_NAVIGATOR
  td->fNsnext++;
  output.AddTrack(track);  
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeomQueryHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector geometry length computation. The tracks are moved into the output basket.
  
// We should make sure that track->fSafety < track->fPstep for these tracks
  constexpr auto kVecSize =  vecCore::VectorSize<vecgeom::VectorBackend::Real_v>();
  TrackVec_t &tracks = input.Tracks();
  const size_t ntr = tracks.size();
// #define DEBUG_VECTOR_DOIT
#ifdef USE_VECGEOM_NAVIGATOR
  // Copy relevant track fields to geometry SOA and process vectorized.
  GeantTrackGeo_v &track_geo = *td->fGeoTrack;
  track_geo.Clear();
  size_t i = 0;
  // Process first the tracks that start from a boundary
  for (size_t itr = 0; itr < ntr; ++itr) {
    GeantTrack *track = tracks[itr];
    // Mark tracks to be processed
    track->fPending = false;
    if (track->fBoundary) {
      td->fPathV[i] = track->Path();
      track_geo.AddTrack(*track, itr);
      track->fPending = true;
      i++;
    }
  }
  if (i >= kVecSize) {
    // Process these tracks in vector mode
    VectorNavInterface::NavFindNextBoundary(i, track_geo.fPstepV,
               track_geo.fXposV, track_geo.fYposV, track_geo.fZposV,
               track_geo.fXdirV, track_geo.fYdirV, track_geo.fZdirV,
               (const VolumePath_t **)td->fPathV,
               track_geo.fSnextV, track_geo.fSafetyV, track_geo.fCompSafetyV);
    for (size_t itr = 0; itr < i; ++itr) {
      // Find the original track
      GeantTrack *track = tracks[track_geo.fIdV[itr]];
#ifdef DEBUG_VECTOR_DOIT
      ScalarNavInterfaceVGM::NavFindNextBoundary(*track);
      if ( vecCore::math::Abs(track->fSnext - track_geo.fSnextV[itr]) > 1.e-10) {
        printf("  track %d (ind=%ld): scalar (snext=%16.12f)  vector (snext=%16.12f)\n",
                track->fParticle, itr, track->fSnext, track_geo.fSnextV[itr]);
      }
#endif
      track->fSnext = track_geo.fSnextV[itr];
      track->fSafety = 0.;
      track->fBoundary = !track_geo.fCompSafetyV[itr];
    }
    track_geo.Clear();
    i = 0;
  }
    
  // Process the tracks having pstep > safety
  for (size_t itr = 0; itr < ntr; ++itr) {
    GeantTrack *track = tracks[itr];
    if (track->fPending || track->fBoundary) continue;
    // Skip tracks with big enough safety
    if (track->fSafety > track->fPstep) {
      track->fSnext = track->fPstep;
      track->fBoundary = false;
      continue;
    }
    td->fPathV[i] = track->Path();
    track_geo.AddTrack(*track, itr);
    track->fPending = true; // probably not needed
    i++;
  }
  if (i >= kVecSize) {  
    VectorNavInterface::NavFindNextBoundary(i, track_geo.fPstepV,
               track_geo.fXposV, track_geo.fYposV, track_geo.fZposV,
               track_geo.fXdirV, track_geo.fYdirV, track_geo.fZdirV,
               (const VolumePath_t **)td->fPathV,
               track_geo.fSnextV, track_geo.fSafetyV, track_geo.fCompSafetyV);  
    for (size_t itr = 0; itr < i; ++itr) {
      // Find the original track
      GeantTrack *track = tracks[track_geo.fIdV[itr]];
#ifdef DEBUG_VECTOR_DOIT
      ScalarNavInterfaceVGM::NavFindNextBoundary(*track);
      if ( vecCore::math::Abs(track->fSnext - track_geo.fSnextV[itr]) > 1.e-10 ||
           vecCore::math::Abs(track->fSafety - track_geo.fSafetyV[itr]) > 1.e-10 ) {
        printf("  track %d (ind=%ld): scalar (snext=%16.12f safety=%16.12f)  vector (snext=%16.12f safety=%16.12f)\n",
                track->fParticle, itr, track->fSnext, track->fSafety, track_geo.fSnextV[itr], track_geo.fSafetyV[itr]);
      }
#endif
      track->fSnext = track_geo.fSnextV[itr];
      track->fSafety = track_geo.fSafetyV[itr];
      track->fBoundary = !track_geo.fCompSafetyV[itr];
    }
  } else {
    for (size_t itr = 0; itr < i; ++itr) {
      GeantTrack *track = tracks[track_geo.fIdV[itr]];
      ScalarNavInterfaceVGM::NavFindNextBoundary(*track);
    }
  }
    
  td->fNsnext += ntr;
  // Copy tracks to output
#ifndef VECCORE_CUDA
  std::move(tracks.begin(), tracks.end(), std::back_inserter(output.Tracks()));
#else
  for (auto track : tracks) output.AddTrack(track);
#endif
#else
// ROOT geometry. Fall back to scalar implementation
  (void)tracks;
  Handler::DoIt(input, output, td);
#endif
}

} // GEANT_IMPL_NAMESPACE
} // Geant
