#include "Geant/Handler.h"
#include "Geant/TaskData.h"
#include "GeantNuma.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler::Handler(int threshold, Propagator *propagator) : fPropagator(propagator)
{
  // Handler constructor. The handler needs to be manually activated to actually
  // allocate the basketizer.
  fBcap = propagator->fConfig->fMaxPerBasket;
  fThreshold.store(threshold);
  // Make sure the threshold is a power of 2
  assert((threshold & (threshold - 1)) == 0 && "Handler threshold must be power of 2");
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fLock.clear();
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler::~Handler()
{
#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
  if (GetNode() < 0)
    delete fBasketizer;
  else
    NumaUtils::NumaAlignedFree(fBasketizer);
#else
  delete fBasketizer;
#endif
  fBasketizer = nullptr;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fLock.clear();
#endif
}

VECCORE_ATT_HOST_DEVICE
Handler::Handler(const Handler &other)
{
  fThreadLocal = other.fThreadLocal;
  fActive = other.fActive;
  fMayBasketize = other.fMayBasketize;
  fId = other.fId;
  fBcap = other.fBcap;
  fThreshold.store(other.fThreshold.load());
  // Do not copy state counters, nor basketizer
  fPropagator = other.fPropagator;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fLock.clear();
#endif
}

VECCORE_ATT_HOST_DEVICE
Handler &Handler::operator=(const Handler &other)
{
  if (&other != this) {
    fThreadLocal = other.fThreadLocal;
    fActive = other.fActive;
    fMayBasketize = other.fMayBasketize;
    fId = other.fId;
    fBcap = other.fBcap;
    fThreshold.store(other.fThreshold.load());
    // Do not copy state counters, nor basketizer
    fPropagator = other.fPropagator;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    fLock.clear();
#endif
  }
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Handler::DoIt(Basket &input, Basket &output, TaskData *td)
{
  // Vector DoIt method implemented as a loop. Overwrite to implement a natively
  // vectorized version.
  for (auto track : input.Tracks())
    DoIt(track, output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Handler::DoItScalar(Basket &input, Basket &output, TaskData *td)
{
  // Basketized DoIt method implemented as a loop. Overwrite to implement a natively
  // vectorized version.
  //  DoIt(input, output, td);
  for (auto track : input.Tracks())
    DoIt(track, output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Handler::ActivateBasketizing(bool flag)
{
  if (fActive == flag) return;
  int basket_size = fThreshold;
  int buffer_size = 1 << 12; // 16 kBytes per handler
  assert(fThreshold < buffer_size);
  // Create basketizer the first time the handler is activated
  if (flag && !fBasketizer) {
#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
    if (GetNode() < 0) {
      fBasketizer = new basketizer_t(buffer_size, basket_size);
    } else {
      int basketizer_size = basketizer_t::SizeofInstance(buffer_size);
      fBasketizer         = basketizer_t::MakeInstanceAt(NumaUtils::NumaAlignedMalloc(basketizer_size, GetNode(), 64),
                                                 buffer_size, basket_size);
    }
#else
    fBasketizer = new basketizer_t(buffer_size, basket_size);
#endif
  }
  fActive = flag;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool Handler::AddTrack(Track *track, Basket &collector, TaskData *)
{
  // Adding a track to the handler assumes that the handler is basketized.
  // The track will be pushed into the basketizer. The calling thread has to
  // provide an empy collector basket which can possibly be filled by the track
  // vector extracted during the operation.

  // Make sure the collector is fit to store the number of tracks required
  collector.Tracks().reserve(fBcap);
  bool extracted = fBasketizer->AddElement(track, collector.Tracks());
  if (extracted) {
    collector.SetThreshold(fThreshold);
  }
  return extracted;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool Handler::Flush(Basket &collector, TaskData *)
{
// Flush if possible remaining tracks from the basketizer into the collector basket.
// NOTE: The operation is not guaranteed to succeed, even if the basketizer
//       contains tracks in case it is 'hot' (e.g. adding tracks or finishing
//       other flush). Flushing is blocking for other flushes!
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  // do not touch if other flushing operation is ongoing
  if (fLock.test_and_set(std::memory_order_acquire)) return false;
#endif
  bool flushed = false;
  while (!flushed && fBasketizer->GetNstored())
    flushed = fBasketizer->Flush(collector.Tracks());
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fLock.clear();
#endif
  return flushed;
}

//______________________________________________________________________________
LocalHandler::LocalHandler(Handler *handler) : Handler(*handler), fHandler(handler)
{
  fThreshold.store(fHandler->GetThreshold());
  fThr = fHandler->GetThreshold();
  fLocalBasket.reserve(fThreshold);
  // Activate basketizing for local handlers that may basketize
  fActive = fMayBasketize;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool LocalHandler::AddTrack(Track *track, Basket &collector, TaskData *)
{
  fLocalBasket.push_back(track);
  if (fLocalBasket.size() == fThr) {
    // Copy tracks in the collector
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    std::copy(fLocalBasket.begin(), fLocalBasket.end(),
              std::back_inserter(collector.Tracks()));
#else
    auto &insertee(td->fStageBuffers[fFollowUpStage]->Tracks());
    for (auto track : fLocalBasket) {
      // If a follow-up stage is declared, this overrides any follow-up set by handlers
      collector.Tracks().push_back(track);
    }
#endif
    fLocalBasket.clear();
    return true;
  }
  return false;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool LocalHandler::Flush(Basket &collector, TaskData *)
{
  bool flush = fLocalBasket.size() > 0;
  if (flush) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    std::copy(fLocalBasket.begin(), fLocalBasket.end(),
              std::back_inserter(collector.Tracks()));
#else
    auto &insertee(td->fStageBuffers[fFollowUpStage]->Tracks());
    for (auto track : fLocalBasket) {
      // If a follow-up stage is declared, this overrides any follow-up set by handlers
      collector.Tracks().push_back(track);
    }
#endif
    fLocalBasket.clear();    
  }
  return flush;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
