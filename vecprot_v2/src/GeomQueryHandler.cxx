#include "GeomQueryHandler.h"
#include "Basket.h"
#include "GeantTaskData.h"
#include "GeantTrackGeo.h"
#include "GeantRunManager.h"
//#include "TransportManager.h"
#include "GeantNuma.h"

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
    if (GetNode() < 0) {
      fBasketizer = new basketizer_t(buffer_size, basket_size);
    } else {
      int basketizer_size = basketizer_t::SizeofInstance(buffer_size);
      fBasketizer = basketizer_t::MakeInstanceAt(
        NumaAlignedMalloc(basketizer_size, GetNode(), 64), buffer_size, basket_size);
    }
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
  // Select follow-up stage
  track->SetStage(ESimulationStage::kPropagationStage);
  output.AddTrack(track);  
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeomQueryHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector geometry length computation. The tracks are moved into the output basket.
  
// We should make sure that track->fSafety < track->fPstep for these tracks
  TrackVec_t &tracks = input.Tracks();
/*
#ifdef USE_VECGEOM_NAVIGATOR
  // Copy relevant track fields to geometry SOA and process vectorized.
  GeantTrackGeo_v &track_geo = *td->fGeoTrack;
  track_geo.Clear();
  int i = 0;
  for (auto track : tracks) {
    track_geo.AddTrack(*track);
    td->fPathV[i] = track->fPath;
    td->fNextpathV[i] = track->fNextpath;
    i++;
  }
    
  // The vectorized SOA call
  VectorNavInterface::NavFindNextBoundary(tracks.size(), track_geo.fPstepV,
               track_geo.fXposV, track_geo.fYposV, track_geo.fZposV,
               track_geo.fXdirV, track_geo.fYdirV, track_geo.fZdirV,
               (const VolumePath_t **)td->fPathV,
               track_geo.fSnextV, track_geo.fSafetyV, track_geo.fBoundaryV);
  
  // Update original tracks
  track_geo.UpdateOriginalTracks();
  // Count this as a single geometry call
  td->fNsnext += 1;
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
*/
  // For the moment just loop and call scalar DoIt
  for (auto track : tracks) {
    DoIt(track, output, td);
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
