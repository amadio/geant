#include "FeederTask.h"

#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif


FeederTask::FeederTask (Geant::GeantTaskData *td): fTd(td) { }

FeederTask::~FeederTask () { }

tbb::task* FeederTask::execute ()
{
  GeantPropagator *propagator = GeantPropagator::Instance();
  int nbaskets = 0;
  if (!propagator->fLastEvent) {
    nbaskets = propagator->ImportTracks(propagator->fNevents, 0, 0, fTd);
    propagator->fLastEvent = propagator->fNevents;
    propagator->ReleaseLock();
    tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
    TransportTask & transportTask = *new(cont.allocate_child()) TransportTask( fTd );
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
      printf("= [task id %d] digitizing event %d with %d tracks pri=%d \n", fTd->fTid, evt->GetEvent(), ntracks, propagator->fPriorityEvents.load());
      //  propagator->fApplication->Digitize(evt->GetEvent());
      propagator->fDoneEvents->SetBitNumber(evt->GetEvent());
      if (propagator->fLastEvent < propagator->fNtotal) {
        printf("=> Importing event %d\n", propagator->fLastEvent);
        nbaskets += propagator->ImportTracks(1, propagator->fLastEvent, islot, fTd);
        propagator->fLastEvent++;
      }
    }
  }

  //printf("=== Feeder task  found %d baskets ===\n", nbaskets);
  // spawn transport task
  propagator->ReleaseLock();
  tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
  TransportTask & transportTask = *new(cont.allocate_child()) TransportTask( fTd );
  return & transportTask;
}
