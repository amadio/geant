#include "Selector.h"
#include "GeantNuma.h"
#include "GeantTaskData.h"
#include "GeantPropagator.h"
#include "GeantRunManager.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Selector::Selector(int threshold, GeantPropagator *propagator,
               int node, int index, Volume_t *vol)
//#ifdef USE_ROOT
//  : TGeoExtension(), fVolume(vol), fIndex(index), fNode(node), fThreshold(threshold) {
//#else
  : fVolume(vol), fIndex(index), fNode(node), fThreshold(threshold) {
//#endif
  // Selector constructor. The selector needs to be manually activated to actually
  // allocate the basketizer.
  fBcap = propagator->fConfig->fMaxPerBasket;
  // Make sure the threshold is a power of 2
  assert(threshold & (threshold - 1) == 0);
  // Connect to logical volume if provided
  if (vol) ConnectToVolume();
  fLock.clear();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Selector::~Selector()
{
  if (fNode < 0) delete fBasketizer;
  else NumaAlignedFree(fBasketizer);
  fBasketizer = nullptr;
  DisconnectVolume();
  fLock.clear();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Selector::ConnectToVolume()
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
void Selector::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector DoIt method implemented as a loop. Overwrite to implement a natively
// vectorized version.
  for (auto track: input.Tracks())
    DoIt(track, output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Selector::DisconnectVolume()
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
void Selector::ActivateBasketizing(bool flag)
{
  if (fActive == flag) return;
  fActive = flag;
  int basket_size = fThreshold;
  int buffer_size = 1 << 12; // 16 kBytes per selector
  // Create basketizer the first time the selector is activated
  if (fActive && ! fBasketizer) {
    // Set a 'compromise' size for the basketizer buffer for all geometry selectors
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
bool Selector::AddTrack(GeantTrack *track, Basket &collector)
{
// Adding a track to the selector assumes that the selector is basketized.
// The track will be pushed into the basketizer. The calling thread has to 
// provide an empy collector basket which can possibly be filled by the track 
// vector extracted during the operation.
  
  // Make sure the collector is fit to store the number of tracks required
  collector.Tracks().reserve(fBcap);
  bool extracted = fBasketizer->AddElement(track, collector.Tracks());
  if (extracted)
    collector.SetThreshold(fThreshold);
  return extracted;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool Selector::Flush(Basket &collector)
{
// Flush if possible remaining tracks from the basketizer into the collector.
// NOTE: The operation is not guaranteed to succeed, even if the basketizer
//       contains tracks in case it is 'hot' (e.g. adding tracks or finishing
//       other flush). Flushing is blocking for other flushes!
  if (fLock.test_and_set(std::memory_order_acquire)) return false;
  bool flushed = fBasketizer->Flush(collector.Tracks());
  fLock.clear();
  return flushed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Basket *Selector::GetFreeBasket(GeantTaskData *td)
{
  Basket *next = td->GetFreeBasket();
  next->Tracks().reserve(fBcap);
  next->SetThreshold(fThreshold);
  return next;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
