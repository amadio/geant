#include "Filter.h"
#include "GeantNuma.h"
#include "GeantTaskData.h"
#include "GeantPropagator.h"
#include "GeantRunManager.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Filter::Filter(int threshold, GeantPropagator *propagator,
               Basket::ELocality locality, int node, int index, Volume_t *vol)
//#ifdef USE_ROOT
//  : TGeoExtension(), fLocality(locality), fVolume(vol), fIndex(index), fNode(node), fThreshold(threshold) {
//#else
  : fLocality(locality), fVolume(vol), fIndex(index), fNode(node), fThreshold(threshold) {
//#endif
  // Filter constructor. The filter needs to be manually activated to actually
  // allocate the basketizer.
  fBcap = propagator->fConfig->fMaxPerBasket;
  // Make sure the threshold is a power of 2
  assert(threshold & (threshold - 1) == 0);
  // Connect to logical volume if provided
  if (vol) ConnectToVolume();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Filter::~Filter()
{
  if (fNode < 0) delete fBasketizer;
  else NumaAlignedFree(fBasketizer);
  fBasketizer = nullptr;
  DisconnectVolume();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Filter::ConnectToVolume()
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
void Filter::DisconnectVolume()
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
void Filter::SetActive(bool flag)
{
  if (fActive == flag) return;
  fActive = flag;
  int basket_size = fThreshold;
  int buffer_size = 1 << 12; // 16 kBytes per filter
  // Create basketizer the first time the filter is activated
  if (fActive && ! fBasketizer) {
    // Set a 'compromise' size for the basketizer buffer for all geometry filters
    if (fVolume) {
      assert(fThreshold < 512);
      buffer_size = 1 << 10; // This makes ~16 MBytes of basketizer buffers for 4K volumes
    }
    assert(fThreshold < buffer_size);
    if (fNode < 0) {
      fBasketizer = new basketizer_t(buffer_size, basket_size);
    } else {
      int basketizer_size = basketizer_t::SizeofInstance(buffer_size);
      fBasketizer = basketizer_t::MakeInstanceAt(
        NumaAlignedMalloc(basketizer_size, fNode, 64), buffer_size, basket_size);
    }
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool Filter::AddTrack(GeantTrack *track, Basket &collector)
{
// Adding a track to the filter will push the track into the basketizer. The 
// calling thread has to provide an empy collector basket which can possibly be
// filled by the track vector extracted during the operation.
  return ( fBasketizer->AddElement(track, collector.Tracks()) );
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool Filter::Flush(Basket &collector)
{
// Flush if possible remaining tracks from the basketizer into the collector.
// NOTE: The operation is not guaranteed to succeed, even if the basketizer
//       contains tracks in case it is 'hot' (e.g. adding tracks or finishing
//       other flush). Flushing is blocking for other flushes!
  return ( fBasketizer->Flush(collector.Tracks()) );
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Basket *Filter::GetFreeBasket(GeantTaskData *td)
{
  Basket *next = td->GetFreeBasket();
  next->SetLocality(fLocality);
  next->Tracks().reserve(fBcap);
  next->SetThreshold(fThreshold);
  return next;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
