#include "Geant/Error.h"

#include "FeederTask.h"
#include "ThreadData.h"
#include "TransportTask.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "GeantEvent.h"
#include "GeantTaskData.h"


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
  if (wm->IsStopped()){
    return NULL;
  }
  GeantPropagator *propagator = GeantPropagator::Instance();
  int tid = wm->ThreadId();
  if(tid>=propagator->fNthreads){
    return NULL;
  }
  Geant::GeantTaskData *td = propagator->fThreadData[tid];
  Geant::Print("","=== Feeder task %d (%d) created ===", tid, td->fTid);


  ThreadData *threadData = ThreadData::Instance(propagator->fNthreads);


  int nbaskets = 0;
  if (!propagator->fLastEvent) {
    nbaskets = propagator->ImportTracks(propagator->fNevents, 0, 0, td);
    propagator->fLastEvent = propagator->fNevents;
    *fNbaskets = nbaskets;
    propagator->ReleaseLock();
    tbb::task::set_ref_count(2);
    TransportTask & transportTask = *new(tbb::task::allocate_child()) TransportTask( nbaskets );
    return & transportTask;
  }
  // Check and mark finished events
  for (int islot = 0; islot < propagator->fNevents; islot++) {
    //Geant::Print("","=== heck and mark finished events %d ===", nbaskets);
    GeantEvent *evt = propagator->fEvents[islot];
    if (propagator->fDoneEvents->TestBitNumber(evt->GetEvent()))
      continue;
    if (evt->Transported()) {
      propagator->fPriorityEvents--;
      evt->Print();
      // Digitizer (todo)
      int ntracks = propagator->fNtracks[islot];
      Printf("= digitizing event %d with %d tracks pri=%d", evt->GetEvent(), ntracks, propagator->fPriorityEvents.load());
      //  propagator->fApplication->Digitize(evt->GetEvent());
      propagator->fDoneEvents->SetBitNumber(evt->GetEvent());
      if (propagator->fLastEvent < propagator->fNtotal) {
        Printf("=> Importing event %d", propagator->fLastEvent);
        nbaskets += propagator->ImportTracks(1, propagator->fLastEvent, islot, td);
        propagator->fLastEvent++;
      }
    }
  }


  *fNbaskets = nbaskets;
  Geant::Print("","=== Feeder task  found %d baskets ===", nbaskets);
  // spawn transport task
  propagator->ReleaseLock();
  tbb::task::set_ref_count(2);
  TransportTask & transportTask = *new(tbb::task::allocate_child()) TransportTask( nbaskets );
  return & transportTask;
}
