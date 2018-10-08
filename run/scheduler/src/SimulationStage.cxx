#include "Geant/SimulationStage.h"

#include "Geant/TaskData.h"
#include "Geant/Propagator.h"
#include "Geant/RunManager.h"
#include "StackLikeBuffer.h"
#include "TrackStat.h"
#include "Geant/BasketCounters.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::SimulationStage(ESimulationStage type, Propagator *prop) : fType(type), fPropagator(prop)
{
  fCheckCountdown = 0;
  fFireFlushRatio = prop->fConfig->fFireFlushRatio;
  fId             = prop->RegisterStage(this);
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fCheckLock.clear();
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::~SimulationStage()
{
  for (int i = 0; i < GetNhandlers(); ++i)
    delete fHandlers[i];
}

VECCORE_ATT_HOST_DEVICE
SimulationStage::SimulationStage(const SimulationStage &other)
{
  fType             = other.fType;
  fPropagator       = other.fPropagator;
  fFireFlushRatio   = other.fFireFlushRatio;
  fId               = other.fId;
  fUserActionsStage = other.fUserActionsStage;
  fFollowUpStage    = other.fFollowUpStage;
  fThrBasketCheck   = other.fThrBasketCheck;
  fNstaged          = other.fNstaged;
  fNbasketized      = other.fNbasketized;
  fUniqueFollowUp   = other.fUniqueFollowUp;
  fEndStage         = other.fEndStage;
  fBasketized       = other.fBasketized;
  for (auto handler : other.fHandlers)
    fHandlers.push_back(handler);
}

VECCORE_ATT_HOST_DEVICE
SimulationStage &SimulationStage::operator=(const SimulationStage &other)
{
  if (&other != this) {
    fType             = other.fType;
    fPropagator       = other.fPropagator;
    fFireFlushRatio   = other.fFireFlushRatio;
    fId               = other.fId;
    fUserActionsStage = other.fUserActionsStage;
    fFollowUpStage    = other.fFollowUpStage;
    fThrBasketCheck   = other.fThrBasketCheck;
    fNstaged          = other.fNstaged;
    fNbasketized      = other.fNbasketized;
    fUniqueFollowUp   = other.fUniqueFollowUp;
    fEndStage         = other.fEndStage;
    fBasketized       = other.fBasketized;
    for (auto handler : other.fHandlers)
      fHandlers.push_back(handler);
  }
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool SimulationStage::HasLocalHandlers() const
{
  for (auto handler : fHandlers) {
    if (handler->IsLocal())
      return true;
  }
  return false;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SimulationStage::ReplaceLocalHandlers()
{
  for (size_t i = 0; i < fHandlers.size(); ++i) {
    if (fHandlers[i]->IsLocal())
      fHandlers[i] = new LocalHandler(fHandlers[i]);
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SimulationStage::DeleteLocalHandlers()
{
  for (size_t i = 0; i < fHandlers.size(); ++i) {
    if (fHandlers[i]->IsLocal()) {
      auto handler = ((LocalHandler*)fHandlers[i])->GetHandler();
      delete fHandlers[i];
      fHandlers[i] = handler;
      fHandlers.clear();
    }
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::CheckBasketizers(TaskData *td, size_t flush_threshold)
{
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  // do not touch if other checking operation is ongoing
  if (fCheckLock.test_and_set(std::memory_order_acquire)) return false;
#endif
  fThrBasketCheck *= 3.;
  fCheckCountdown = fThrBasketCheck;
  int nactivated  = 0;
  size_t nactive  = 0;
  Basket &output  = *td->fShuttleBasket;
  output.Clear();
  size_t nflushedtot = 0;
  size_t nfiredtot   = 0;
  for (int i = 0; i < GetNhandlers(); ++i) {
    if (!fHandlers[i]->MayBasketize()) continue;
    size_t nflushed = fPropagator->fRunMgr->GetNflushed(fId, i);
    size_t nfired   = fPropagator->fRunMgr->GetNflushed(fId, i);
    nflushedtot += nflushed;
    nfiredtot += nfired;
    if (fHandlers[i]->IsActive()) {
      // Now we check if the firing rate compared to flushing is too low.
      // We could first try to reduce the basket size.
      nactive++;
      if (nfired < nflushed * flush_threshold) {
        fHandlers[i]->ActivateBasketizing(false);
        FlushHandler(i, td, output);
        nactivated--;
        nactive--;
      }
    } else {
      // We start basketizing if firing rate is above threshold
      if (nfired > nflushed * flush_threshold) {
        fHandlers[i]->ActivateBasketizing(true);
        nactivated++;
        nactive++;
      }
    }
  }
  float nbasketized = td->fCounters[fId]->fNvector;
  float ntotal      = nbasketized + td->fCounters[fId]->fNscalar;
  Printf("Stage %20s: basketized %d %% | nscalar = %ld  nvector = %ld  nfired = %ld nflushed = %ld", GetName(),
         int(100 * nbasketized / ntotal), size_t(ntotal - nbasketized), size_t(nbasketized), nfiredtot, nflushedtot);

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fCheckLock.clear();
#endif
  CopyToFollowUps(output, td);
  if (nactivated > 0)
    Printf("--- activated %d basketizers for stage %s : now %ld/%ld active\n", nactivated, GetName(), nactive,
           fHandlers.size());
  else if (nactivated < 0)
    Printf("--- stopped   %d basketizers for stage %s : now %ld/%ld active\n", -nactivated, GetName(), nactive,
           fHandlers.size());
  return nactivated;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::FlushHandler(int i, TaskData *td, Basket &output)
{
  if (!fHandlers[i]->HasTracks()) return 0;
  Basket &bvector = *td->fBvector;
  bvector.Clear();
  fHandlers[i]->Flush(bvector, td);
  td->fCounters[fId]->fFlushed[i]++;
  int nflushed = bvector.size();
  if (nflushed >= fPropagator->fConfig->fNvecThreshold) {
    td->fCounters[fId]->fNvector += bvector.size();
    fHandlers[i]->DoIt(bvector, output, td);
  } else {
    td->fCounters[fId]->fNscalar += bvector.size();
    for (auto track : bvector.Tracks())
      fHandlers[i]->DoIt(track, output, td);
  }
  return nflushed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::FlushAndProcess(TaskData *td)
{
  // Flush all active handlers in the stage, executing their scalar DoIt methods.
  // Flushing the handlers is opportunistic, as a handler is flushed by the first
  // thread to come. Subsequent threads will just notice that the handler is busy
  // and continue to other handlers, then to the next stage.

  Basket &input = *td->fStageBuffers[fId];
  if (fBasketized) {
    int check_countdown = fCheckCountdown.fetch_sub(input.Tracks().size()) - input.Tracks().size();
    if (check_countdown <= 0) CheckBasketizers(td, fFireFlushRatio);
  }
  Basket &bvector = *td->fBvector;
  bvector.Clear();
  Basket &output = *td->fShuttleBasket;
  output.Clear();
  // Loop tracks in the input basket and select the appropriate handler
  for (auto track : input.Tracks()) {
    // Default next stage is the follow-up
    track->SetStage(fFollowUpStage);
    Handler *handler = Select(track, td);
    // Execute in scalar mode the handler action
    if (handler) {
      td->fCounters[fId]->Increment(handler->GetId(), handler->GetThreshold());
      td->fCounters[fId]->fNscalar += size_t(fBasketized & handler->MayBasketize());
      handler->DoIt(track, output, td);
    } else
      output.AddTrack(track);
  }
  // The stage buffer needs to be cleared
  input.Clear();

  // Loop active handlers and flush them into btodo basket
  int nflush = 0;
  for (int i = 0; i < GetNhandlers(); ++i) {
    int nbasket = 0;
    if (fHandlers[i]->IsActive() && (nbasket = fHandlers[i]->Flush(bvector, td))) {
      // btodo has some content, invoke DoIt
      td->fCounters[fId]->fFlushed[i]++;
      if (bvector.size() >= (size_t)fPropagator->fConfig->fNvecThreshold) {
        td->fCounters[fId]->fNvector += bvector.size();
        fHandlers[i]->DoIt(bvector, output, td);
      } else {
        td->fCounters[fId]->fNscalar += bvector.size();
        for (auto track : bvector.Tracks())
          fHandlers[i]->DoIt(track, output, td);
      }
      bvector.Clear();
    }
    nflush += nbasket;
    // if (nflush > fPropagator->fConfig->fNvecThreshold) break;
  }

  return CopyToFollowUps(output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::Process(TaskData *td)
{
  // Processing is concurrent for all tasks/threads serving the same propagator.
  // The input basket is the task data-specific container for the caller thread
  // corresponding to this stage. The DoIt method is executed for the
  // handlers present in the stage, in scalar mode for non-active ones and using the
  // vector interface for the active ones. The resulting tracks after the stage
  // execution are transported via a thread-specific shuttle basket into the inputs
  // of the follow-up stages set by the handler DoIt method. The method returns
  // the number of tracks pushed to the output.

  constexpr int kMinGen = 0; // minimum generation
  assert(fFollowUpStage >= 0);
  Basket &input = *td->fStageBuffers[fId];
  int ninput    = input.Tracks().size();
  if (fBasketized) {
    int check_countdown = fCheckCountdown.fetch_sub(ninput) - ninput;
    if (check_countdown <= 0) CheckBasketizers(td, fFireFlushRatio);
  }
  Basket &bvector = *td->fBvector;
  bvector.Clear();
  Basket &output = *td->fShuttleBasket;
  output.Clear();

#ifdef GEANT_DEBUG
// static long counter = 0;
// counter++;
#endif

  // Loop tracks in the input basket and select the appropriate handler
  for (auto track : input.Tracks()) {
// Default next stage is the follow-up
#ifdef GEANT_DEBUG
//    assert(!input.HasTrackMany(track));
//    assert(!output.HasTrack(track));
#endif
    track->SetStage(fFollowUpStage);
    Handler *handler = Select(track, td);
    // If no handler is selected the track does not perform this stage
    if (!handler) {
      output.AddTrack(track);
      continue;
    }

    if (!handler->IsActive() || track->GetGeneration() < kMinGen) {
      if (fBasketized) {
        td->fCounters[fId]->Increment(handler->GetId(), handler->GetThreshold());
        // Account only scalar calls to basketizable handlers
        if (handler->MayBasketize()) td->fCounters[fId]->fNscalar++;
      }
      // Scalar DoIt.
      handler->DoIt(track, output, td);
    } else {
      // Add the track to the handler, which may extract a full vector.
      if (handler->AddTrack(track, bvector, td)) {
        // Vector DoIt
        td->fCounters[fId]->fNvector += bvector.size();
        td->fCounters[fId]->fFired[handler->GetId()]++;
        handler->DoIt(bvector, output, td);
        bvector.Clear();
      }
    }
  }
  // The stage buffer needs to be cleared
  input.Clear();
  // The track and its eventual progenies should be now copied to the output
  return CopyToFollowUps(output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::CopyToFollowUps(Basket &output, TaskData *td)
{
  // Copy tracks from output basket to follow-up stages. Output needs to be cleared.
  int ntracks = output.size();

  if (fEndStage) {
    // The last stage copies all tracks into the stage buffer
    td->fStackBuffer->AddTracks(output.Tracks());
    return ntracks;
  }
  // Copy output tracks to the follow-up stages
  if (fUniqueFollowUp) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    std::copy(output.Tracks().begin(), output.Tracks().end(),
              std::back_inserter(td->fStageBuffers[fFollowUpStage]->Tracks()));
#else
    auto &insertee(td->fStageBuffers[fFollowUpStage]->Tracks());
    for (auto track : output.Tracks()) {
      // If a follow-up stage is declared, this overrides any follow-up set by handlers
      insertee.push_back(track);
    }
#endif
  } else {
    for (auto track : output.Tracks()) {
      assert(track->GetStage() != fId); // no stage feeds itself
      td->fStageBuffers[track->GetStage()]->AddTrack(track);
    }
  }
  return ntracks;
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
