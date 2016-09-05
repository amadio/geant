#include "FeederTask.h"

#include "GeantConfig.h"
#include "GeantRunManager.h"

#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif


FeederTask::FeederTask (Geant::GeantTaskData *td): fTd(td) { }

FeederTask::~FeederTask () { }

tbb::task* FeederTask::execute ()
{
  GeantPropagator *propagator = fTd->fPropagator;
  GeantRunManager *runmgr = propagator->fRunMgr;
  int nbaskets = 0;
  if (!propagator->fConfig->fLastEvent) {
    nbaskets = runmgr->ImportTracks(propagator->fConfig->fNbuff, 0, 0, fTd);
    propagator->fConfig->fLastEvent = propagator->fConfig->fNbuff;
    runmgr->ReleaseLock();
    tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
    TransportTask & transportTask = *new(cont.allocate_child()) TransportTask( fTd );
    return & transportTask;
  }
  // Check and mark finished events
  for (int islot = 0; islot < propagator->fConfig->fNbuff; islot++) {
    //Geant::Print("","=== check and mark finished events %d ===", nbaskets);
    GeantEvent *evt = runmgr->GetEvent(islot);
    if (runmgr->IsDoneEvent(evt->GetEvent()))
      continue;
    if (evt->Transported()) {
      (runmgr->GetPriorityEvents())--;
      evt->Print();
      // Digitizer (todo)
      int ntracks = runmgr->GetNtracks(islot);
      printf("= [task id %d] digitizing event %d with %d tracks pri=%d \n", fTd->fTid, evt->GetEvent(), ntracks, runmgr->GetNpriority());
      //  propagator->fApplication->Digitize(evt->GetEvent());
      runmgr->SetDoneEvent(evt->GetEvent());
      if (propagator->fConfig->fLastEvent < propagator->fConfig->fNtotal) {
        printf("=> Importing event %d\n", propagator->fConfig->fLastEvent);
        nbaskets += runmgr->ImportTracks(1, propagator->fConfig->fLastEvent, islot, fTd);
        propagator->fConfig->fLastEvent++;
      }
    }
  }

  //printf("=== Feeder task  found %d baskets ===\n", nbaskets);
  // spawn transport task
  runmgr->ReleaseLock();
  tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
  TransportTask & transportTask = *new(cont.allocate_child()) TransportTask( fTd );
  return & transportTask;
}
