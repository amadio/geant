#include "SimulationStage.h"

#include "GeantTaskData.h"
#include "GeantPropagator.h"
#include "StackLikeBuffer.h"
#include "TrackStat.h"
#include "BasketCounters.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::SimulationStage(ESimulationStage type, GeantPropagator *prop)
  : fType(type), fPropagator(prop)
{
  fCheckCountdown = 0;
  fFireFlushRatio = prop->fConfig->fFireFlushRatio;
  fId = prop->RegisterStage(this);
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fCheckLock.clear();
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::~SimulationStage()
{
  for (int i=0; i<GetNhandlers(); ++i)
    delete fHandlers[i];  
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::CheckBasketizers(GeantTaskData *td, size_t flush_threshold)
{
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    // do not touch if other checking operation is ongoing
  if (fCheckLock.test_and_set(std::memory_order_acquire)) return false;
#endif
  fThrBasketCheck *= 2.;
  fCheckCountdown = fThrBasketCheck;
  int nactivated = 0;
  size_t nactive = 0;
  Basket &output = *td->fShuttleBasket;
  output.Clear();
  size_t nflushedtot = 0;
  size_t nfiredtot = 0;
  for (int i=0; i<GetNhandlers(); ++i) {
    if (!fHandlers[i]->MayBasketize()) continue;
    size_t nflushed = fHandlers[i]->GetNflushed();
    size_t nfired = fHandlers[i]->GetNfired();
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
      }
    } else {
      // We start basketizing if firing rate is above threshold
      if (nfired > nflushed * flush_threshold) {
        fHandlers[i]->ActivateBasketizing(true);
        nactivated++;
      }
    }
  }
  float nbasketized = td->fCounters[fId]->fNvector;
  float ntotal = nbasketized + td->fCounters[fId]->fNscalar;
  Printf("Stage %20s: basketized %d %% | nscalar = %ld  nvector = %ld  nfired = %ld nflushed = %ld",
        GetName(), int(100*nbasketized/ntotal),
        size_t(ntotal-nbasketized), size_t(nbasketized), nfiredtot, nflushedtot);
  
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fCheckLock.clear();
#endif
  CopyToFollowUps(output, td);
  if (nactivated > 0)      Printf("--- activated %d basketizers for stage %s\n", nactivated, GetName());
  else if (nactivated < 0) Printf("--- stopped %d basketizers for stage %s\n", -nactivated, GetName());
  printf("%ld/%ld active basketizers\n", nactive, fHandlers.size());
  return nactivated;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::FlushHandler(int i, GeantTaskData *td, Basket &output)
{
  if (!fHandlers[i]->HasTracks()) return 0;
  Basket &bvector = *td->fBvector;
  bvector.Clear();
  fHandlers[i]->Flush(bvector);
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
int SimulationStage::FlushAndProcess(GeantTaskData *td)
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
      td->fCounters[fId]->fNscalar++;
      handler->DoIt(track, output, td);
    }
    else output.AddTrack(track);
  }
  // The stage buffer needs to be cleared
  input.Clear();
  
  // Loop active handlers and flush them into btodo basket
  for (int i=0; i < GetNhandlers(); ++i) {
    if (fHandlers[i]->IsActive() && fHandlers[i]->Flush(bvector)) {
      // btodo has some content, invoke DoIt
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
  }
    
  return CopyToFollowUps(output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::Process(GeantTaskData *td)
{
// Processing is concurrent for all tasks/threads serving the same propagator.
// The input basket is the task data-specific container for the caller thread
// corresponding to this stage. The DoIt method is executed for the 
// handlers present in the stage, in scalar mode for non-active ones and using the
// vector interface for the active ones. The resulting tracks after the stage
// execution are transported via a thread-specific shuttle basket into the inputs
// of the follow-up stages set by the handler DoIt method. The method returns
// the number of tracks pushed to the output.

  assert(fFollowUpStage >= 0);
  Basket &input = *td->fStageBuffers[fId];
  int ninput = input.Tracks().size();
  if (fBasketized) {
    int check_countdown = fCheckCountdown.fetch_sub(ninput) - ninput;
    if (check_countdown <= 0) CheckBasketizers(td, fFireFlushRatio);
  }
  Basket &bvector = *td->fBvector;
  bvector.Clear();
  Basket &output = *td->fShuttleBasket;
  output.Clear();

  #ifdef GEANT_DEBUG
  //static long counter = 0;
  //counter++;
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
    
    if (!handler->IsActive()) {
      // Scalar DoIt.
      // The track and its eventual progenies should be now copied to the output
      if (fBasketized) {
        size_t &counter = td->fCounters[fId]->fCounters[handler->GetId()];
        counter++;
        if ((int)counter == handler->GetThreshold()) {
          counter = 0;
          handler->AddFired();
        }
      }
      td->fCounters[fId]->fNscalar++;
      handler->DoIt(track, output, td);
      #ifdef GEANT_DEBUG
//      assert(!output.HasTrackMany(track));
      #endif
    } else {
      // Add the track to the handler, which may extract a full vector.
      if (handler->AddTrack(track, bvector)) {
      // Vector DoIt
        td->fCounters[fId]->fNvector += bvector.size();
        handler->DoIt(bvector, output, td);
        bvector.Clear();
      }
    }
  }
  // The stage buffer needs to be cleared
  input.Clear();
  return CopyToFollowUps(output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::CopyToFollowUps(Basket &output, GeantTaskData *td)
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
    auto &insertee( td->fStageBuffers[fFollowUpStage]->Tracks() );
    for (auto track : output.Tracks()) {
      // If a follow-up stage is declared, this overrides any follow-up set by handlers
      insertee.push_back(track);
    }
#endif
  } else {
    for (auto track : output.Tracks()) {
      assert(track->GetStage() != fId);         // no stage feeds itself
      td->fStageBuffers[track->GetStage()]->AddTrack(track);
    }
  }
  return ntracks;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
