#include "SimulationStage.h"

#include "Handler.h"
#include "GeantTaskData.h"
#include "GeantPropagator.h"
#include "StackLikeBuffer.h"
#include "TrackStat.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::SimulationStage(ESimulationStage type, GeantPropagator *prop)
  : fType(type), fPropagator(prop)
{
  fId = prop->RegisterStage(this);
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
int SimulationStage::FlushAndProcess(GeantTaskData *td)
{
// Flush all active handlers in the stage, executing their scalar DoIt methods.
// Flushing the handlers is opportunistic, as a handler is flushed by the first
// thread to come. Subsequent threads will just notice that the handler is busy
// and continue to other handlers, then to the next stage.

  Basket &input = *td->fStageBuffers[fId];
  Basket &bvector = *td->fBvector;
  bvector.Clear();
  Basket &output = *td->fShuttleBasket;
  output.Clear();
  int ninput = 0;
  // Loop tracks in the input basket and select the appropriate handler
  for (auto track : input.Tracks()) {
    Handler *handler = Select(track, td);
    // Execute in scalar mode the handler action
    if (handler)
      handler->DoIt(track, output, td);
    else output.AddTrack(track);
  }
  ninput += input.size();
  // The stage buffer needs to be cleared
  input.Clear();
  
  // Loop active handlers and flush them into btodo basket
  for (int i=0; i < GetNhandlers(); ++i) {
    if (fHandlers[i]->IsActive() && fHandlers[i]->Flush(bvector)) {
      // btodo has some content, invoke DoIt
      if (bvector.size() >= fPropagator->fConfig->fNvecThreshold) {
        fHandlers[i]->DoIt(bvector, output, td);
      } else {      
        for (auto track : bvector.Tracks())
          fHandlers[i]->DoIt(track, output, td);
      }
      ninput += bvector.size();
      bvector.Clear();
    }
  }
  
  td->fStat->AddTracks(output.size() - ninput);
  
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

  Basket &input = *td->fStageBuffers[fId];
  Basket &bvector = *td->fBvector;
  bvector.Clear();
  Basket &output = *td->fShuttleBasket;
  output.Clear();
  int ninput = 0;
  // Loop tracks in the input basket and select the appropriate handler
  for (auto track : input.Tracks()) {
    Handler *handler = Select(track, td);
    // If no handler is selected the track does not perform this stage
    if (!handler) {
      output.AddTrack(track);
      ninput++;
      continue;
    }
    
    if (!handler->IsActive()) {
      // Scalar DoIt.
      // The track and its eventual progenies should be now copied to the output
      handler->DoIt(track, output, td);
      ninput++;
    } else {
      // Add the track to the handler, which may extract a full vector.
      if (handler->AddTrack(track, bvector)) {
      // Vector DoIt
        ninput += bvector.size();
        handler->DoIt(bvector, output, td);
        bvector.Clear();
      }
    }
  }
  // The stage buffer needs to be cleared
  input.Clear();
  td->fStat->AddTracks(output.size()-ninput);
  return CopyToFollowUps(output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::CopyToFollowUps(Basket &output, GeantTaskData *td)
{
// Copy tracks from output basket to follow-up stages. Output needs to be cleared.
  int ntracks = output.size();
  
  if (fEndStage) {
    td->fStackBuffer->AddTracks(output.Tracks());
    return ntracks;
  }
  // Copy output tracks to the follow-up stages
  if (fFollowUpStage) {
    for (auto track : output.Tracks()) {
      // If a follow-up stage is declared, this overrides any follow-up set by handlers
      track->fStage = fFollowUpStage;
    }
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
      assert(track->fStage != fId);         // no stage feeds itself
      td->fStageBuffers[track->fStage]->AddTrack(track);
    }
  }
  return ntracks;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
