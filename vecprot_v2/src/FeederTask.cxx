#include "FeederTask.h"

#include <iostream>
#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif

FeederTask::FeederTask (Geant::GeantTaskData *td, int * nbaskets): fTd(td), fNbaskets(nbaskets) { }

FeederTask::~FeederTask () { }

tbb::task* FeederTask::execute ()
{

  GeantPropagator *propagator = GeantPropagator::Instance();

  // Task spawned to inject the next event(s)
  // Only one task at a time
  if (propagator->fFeederLock.test_and_set(std::memory_order_acquire)){
    *fNbaskets = -1;
    return NULL;
  }
  int nbaskets = 0;
  if (!propagator->fLastEvent) {
    nbaskets = propagator->ImportTracks(propagator->fNevents, 0, 0, fTd);
    propagator->fLastEvent = propagator->fNevents;
    propagator->fFeederLock.clear(std::memory_order_release);
    *fNbaskets = nbaskets;
    return NULL;
  }
  // Check and mark finished events
  for (int islot = 0; islot < propagator->fNevents; islot++) {
    GeantEvent *evt = propagator->fEvents[islot];
    if (propagator->fDoneEvents->TestBitNumber(evt->GetEvent()))
      continue;
    if (evt->Transported()) {
      propagator->fPriorityEvents--;
      evt->Print();
      // Digitizer (todo)
      int ntracks = propagator->fNtracks[islot];
      Printf("= digitizing event %d with %d tracks pri=%d", evt->GetEvent(), ntracks, propagator->fPriorityEvents.load());
      //            propagator->fApplication->Digitize(evt->GetEvent());
      propagator->fDoneEvents->SetBitNumber(evt->GetEvent());
      if (propagator->fLastEvent < propagator->fNtotal) {
        Printf("=> Importing event %d", propagator->fLastEvent);
        nbaskets += propagator->ImportTracks(1, propagator->fLastEvent, islot, fTd);
        propagator->fLastEvent++;
      }
    }
  }

  propagator->fFeederLock.clear(std::memory_order_release);
  *fNbaskets = nbaskets;

  return NULL;
}
