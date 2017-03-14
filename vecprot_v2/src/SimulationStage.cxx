#include "SimulationStage.h"

#include "Handler.h"
#include "GeantTaskData.h"
#include "GeantPropagator.h"
#include "StackLikeBuffer.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::SimulationStage(ESimulationStage type, GeantPropagator *prop)
  : fType(type), fPropagator(prop)
{
  CreateHandlers();
  assert((GetNhandlers() > 0) && "Number of handlers for a simulation stage cannot be 0");
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
  Basket &bvector = *td->fBvector;
  Basket &output = *td->fShuttleBasket;
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

  Basket &input = *td->fStageBuffers[fId];
  Basket &bvector = *td->fBvector;
  Basket &output = *td->fShuttleBasket;
// Loop tracks in the input basket and select the appropriate handler
  for (auto track : input.Tracks()) {
    Handler *handler = Select(track, td);
    // If no handler is selected the track does not perform this stage
    if (!handler) {
      output.AddTrack(track);
      continue;
    }
    
    if (!handler->IsActive()) {
      // Scalar DoIt.
      // The track and its eventual progenies should be now copied to the output
      handler->DoIt(track, output, td);
    } else {
      // Add the track to the handler, which may extract a full vector.
      if (handler->AddTrack(track, bvector)) {
      // Vector DoIt
        handler->DoIt(bvector, output, td);
        // The tracks from bvector are now copied into the output
        bvector.Clear();
      }
    }
  }
  input.Clear();
  return CopyToFollowUps(output, td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::CopyToFollowUps(Basket &output, GeantTaskData *td)
{
// Copy tracks from output basket to follow-up stages
  int ntracks = output.size();
  
  if (fEndStage) {
    td->fStackBuffer->AddTracks(output.Tracks());
    return ntracks;
  }
  // Copy output tracks to the follow-up stages
  if (fFollowUpStage) {
    for (auto track : output.Tracks()) {
      // If a follow-up stage is declared, this overrides any follow-up set by handlers
      if (fFollowUpStage) track->fStage = fFollowUpStage;
    }
    std::copy(output.Tracks().begin(), output.Tracks().end(),
              std::back_inserter(td->fStageBuffers[fFollowUpStage]->Tracks()));
  } else {    
    for (auto track : output.Tracks()) {
      assert(track.fStage != fId);         // no stage feeds itself
      td->fStageBuffers[track->fStage]->AddTrack(track);
    }
  }
  output.Clear();
  return ntracks;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
