#include "Selector.h"
#include "GeantNuma.h"
#include "GeantTaskData.h"
#
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Selector::Selector(int threshold, GeantPropagator *propagator)
  : fPropagator(propagator) {
  // Selector constructor. The selector needs to be manually activated to actually
  // allocate the basketizer.
  fBcap = propagator->fConfig->fMaxPerBasket;
  fThreshold.store(threshold);
  // Make sure the threshold is a power of 2
  assert(threshold & (threshold - 1) == 0 && "Selector threshold must be power of 2");
  fLock.clear();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Selector::~Selector()
{
  if (GetNode() < 0) delete fBasketizer;
  else NumaAlignedFree(fBasketizer);
  fBasketizer = nullptr;
  fLock.clear();
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
void Selector::ActivateBasketizing(bool flag)
{
  if (fActive == flag) return;
  fActive = flag;
  int basket_size = fThreshold;
  int buffer_size = 1 << 12; // 16 kBytes per selector
  assert(fThreshold < buffer_size);
  // Create basketizer the first time the selector is activated
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
