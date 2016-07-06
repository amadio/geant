#include "Geant/Error.h"

#include "FeederTask.h"
#include "ThreadData.h"
#include "TransportTask.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "GeantEvent.h"


#include <iostream>
#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif

using namespace Geant;

FeederTask::FeederTask (int *nbaskets): fNbaskets(nbaskets) { }

FeederTask::~FeederTask () { }

tbb::task* FeederTask::execute ()
{
  WorkloadManager *wm = WorkloadManager::Instance();
  int tid = wm->Instance()->ThreadId();
  printf("Feedertask created");
  GeantPropagator *propagator = GeantPropagator::Instance();
  GeantTaskData *td = propagator->fThreadData[tid];
  ThreadData *threadData = ThreadData::Instance(propagator->fNthreads);

  // Task spawned to inject the next event(s)
  // Only one task at a time
  if (propagator->fFeederLock.test_and_set(std::memory_order_acquire)){
    *fNbaskets = -1;
    return NULL;
  }
  int nbaskets = 0;
  if (!propagator->fLastEvent) {
    nbaskets = propagator->ImportTracks(propagator->fNevents, 0, 0, td);
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
        nbaskets += propagator->ImportTracks(1, propagator->fLastEvent, islot, td);
        propagator->fLastEvent++;
      }
    }
  }

  propagator->fFeederLock.clear(std::memory_order_release);
  *fNbaskets = nbaskets;

  // spawn transport task
  TransportTask & transportTask = *new(tbb::task::allocate_child()) TransportTask();
  tbb:task::spawn(transportTask);

  return NULL;
}
