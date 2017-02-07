#include "SimulationStage.h"

#include "Filter.h"
#include "GeantTaskData.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::SimulationStage(EStage type) : fType(type)
{
  fToProcess = new queue_t(1<<16);
  fFilters = new Filter*[fCapacity];
  fNfilters = CreateFilters();
  assert((fNfilters > 0) && "Number of filters for a simulation stage cannot be 0");
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void SimulationStage::AddFilter(Filter *filter)
{
  if (fNfilters == fCapacity) {
    fCapacity *= 2;
    Filter **filters = new Filter*[fCapacity];
    memcpy(filters, fFilters, fNfilters * sizeof(Filter*));
    delete [] fFilters;
    fFilters = filters;
  }
  fFilters[fNfilters++] = filter;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage::~SimulationStage()
{
  for (int i=0; i<fNfilters; ++i)
    delete fFilters[i];
  delete [] fFilters;
  delete fToProcess;
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
  Basket *output = nullptr;
  for (int i=0; i < fNfilters; ++i) {
    if (fFilters[i]->IsActive() && fFilters[i]->Flush(btodo)) {
      if (!output) output = td->GetFreeBasket();
      for (auto track : btodo.Tracks())
        fFilters[i]->DoIt(track, *output, td);
      btodo.Clear();
      int nout = output->GetTracks().size();
      ntracks += nout;
      if (nout > output->GetThreshold()) {
        // Copy now the output to the next stage
        output->SetStage(fNextStage);
        if (!td->fLastBasket)
          td->fLastBasket = output;
        else
          fNextStage->AddBasket(output);
        output = nullptr;
      }
    }
  }
  if (output) {
    if (!output->GetTracks().size())
      td->RecycleBasket(output);
    else {
      output->SetStage(fNextStage);
      if (!td->fLastBasket)
        td->fLastBasket = output;
      else
        fNextStage->AddBasket(output);
    }
  }
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
  Basket *bvector = nullptr;
  Basket *output = td->GetFreeBasket();
// Loop tracks in the basket and select the appropriate filter
  for (auto track : input.Tracks()) {
    Filter *filter = Select(track);
    // If no filter is selected the track does not perform this stage
    if (!filter) output->AddTrack(track);
    if (!filter->IsActive()) {
      // Scalar DoIt
      filter->DoIt(track, *output, td);
      ntracks++;
    } else {
      // Add the track to the filter, which may extract a full vector.
      if (!bvector) bvector = td->GetFreeBasket();
      if (filter->AddTrack(track, *bvector)) {
      // Vector DoIt
        ntracks += bvector->GetTracks().size();
        filter->DoIt(*bvector, *output, td);
      }
    }
  }
  // Recycle used basket
  if (bvector) {
    bvector->Clear();
    td->RecycleBasket(bvector);
  }
  // Mark output for the next stage. Keep the basket in td or push it into
  // next stage queue.
  output->SetStage(fNextStage);
  if (!td->fLastBasket)
    td->fLastBasket = output;
  else
    fNextStage->AddBasket(output);
  return ntracks;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
