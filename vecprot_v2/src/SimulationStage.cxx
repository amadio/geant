#include "SimulationStage.h"

#include "Filter.h"
#include "GeantTaskData.h"
#include "GeantPropagator.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::SimulationStage(ESimulationStage type, GeantPropagator *prop)
  : fType(type), fPropagator(prop)
{
  CreateFilters();
  assert((GetNfilters() > 0) && "Number of filters for a simulation stage cannot be 0");
  fId = prop->RegisterStage(this);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SimulationStage::AddFilter(Filter *filter)
{
  fFilters.push_back(filter);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::~SimulationStage()
{
  for (int i=0; i<GetNfilters(); ++i)
    delete fFilters[i];  
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::FlushAndProcess(Basket &btodo, GeantTaskData *td)
{
// Flush all active filters in the stage, executing their scalar DoIt methods.
// Flushing the filters is opportunistic, as a filter is flushed by the first
// thread to come. Subsequent threads will just notice that the filter is busy
// and continue to other filters, then to the next stage.
  int ntracks = 0;
  btodo.SetStage(this);
  // In case the follow-up stage is unique, the output can be dumped into
  // its buffer.
  Basket &output = (NFollowUps() == 1) ? *td->fStageBuffers[fId]
                                       : *td->GetFreeBasket();
  // Loop active filters and flush them into btodo basket
  for (int i=0; i < GetNfilters(); ++i) {
    if (fFilters[i]->IsActive() && fFilters[i]->Flush(btodo)) {
      // btodo has some content, invoke scalar DoIt
      ntracks += btodo.GetTracks().size();
      for (auto track : btodo.Tracks())
        fFilters[i]->DoIt(track, output, td);
      btodo.Clear();
      // If more follow-up stages, copy to the right buffer
      if (NFollowUps() > 1) {
        for (auto processed : output.Tracks())
          td->fStageBuffers[SelectFollowUp(processed)]->AddTrack(processed);
        output.Clear();
      }
    }
  }
  if (NFollowUps() > 1)
    td->RecycleBasket(&output);
  return ntracks;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SimulationStage::Process(Basket &input, GeantTaskData *td)
{
// A task/thread having the task data td processed a basket for this stage.
// Processing is concurrent for all tasks/threads serving the same propagator.
// The input basket is normally taken from the input queue for the stage, but
// it can be coming from other source. The DoIt method is executed for the 
// filters present in the stage, in scalar mode for non-active ones and using the
// vector interface for the active ones. The resulting tracks after the stage
// execution have to be filled in the output basket extracted from the task data
// pool, which is then either kept to be processed by the same task/thread or
// pushed into the next stage queue.
  int ntracks = 0;
  Basket &bvector = *td->fBvector;
  // In case the follow-up stage is unique, the output can be dumped into
  // its buffer.
  Basket &output = (NFollowUps() == 1) ? *td->fStageBuffers[fId]
                                       : *td->GetFreeBasket();
// Loop tracks in the input basket and select the appropriate filter
  for (auto track : input.Tracks()) {
    Filter *filter = SelectFilter(track);
    // If no filter is selected the track does not perform this stage
    if (!filter) {
      output.AddTrack(track);
      continue;
    }
    
    if (!filter->IsActive()) {
      // Scalar DoIt.
      // The track and its eventual progenies should be now copied to the output
      filter->DoIt(track, output, td);
      ntracks++;
    } else {
      // Add the track to the filter, which may extract a full vector.
      if (filter->AddTrack(track, bvector)) {
      // Vector DoIt
        ntracks += bvector.GetTracks().size();
        filter->DoIt(bvector, output, td);
      }
    }
  }
  bvector.Clear();
  // Mark output for the next stage. Keep the basket in td or push it into
  // next stage queue.
  if (NFollowUps() == 1) return ntracks;

  for (auto track : output.Tracks())
    td->fStageBuffers[SelectFollowUp(track)]->AddTrack(track);

  output.Clear();
  td->RecycleBasket(&output);

  return ntracks;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
