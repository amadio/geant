#include "FlowControllerTask.h"

#include "ThreadData.h"
#include "FeederTask.h"
#include "GeantRunManager.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"

#ifdef USE_ROOT
#include "TThreadMergingFile.h"
#endif

#include "tbb/task_scheduler_init.h"


FlowControllerTask::FlowControllerTask (Geant::GeantTaskData *td, bool starting): fTd(td), fStarting(starting) { }

FlowControllerTask::~FlowControllerTask () { }

tbb::task* FlowControllerTask::execute ()
{
  GeantPropagator *propagator = fTd->fPropagator;
  GeantRunManager *runmgr = propagator->fRunMgr;
  WorkloadManager *wm = propagator->fWMgr;

  //printf("=== %d Flow Controller  ===\n", fTd->fTid);
  if(fStarting){
    while(runmgr->TryLock())
      ;

    tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
    FeederTask & feederTask = *new(cont.allocate_child()) FeederTask( fTd );
    return & feederTask;
  }

  ThreadData *threadData = ThreadData::Instance(runmgr->GetNthreadsTotal());
#ifdef USE_ROOT
  Geant::TThreadMergingFile* file = threadData->fFiles[fTd->fTid];
  bool concurrentWrite = propagator->fConfig->fConcurrentWrite &&
                         propagator->fConfig->fFillTree;
#endif

  if (runmgr->TransportCompleted()) {
    //finish tasks
    //tbb::task::destroy(*tbb::task::parent());
    printf("=== Exit thread %d from Flow Controller  ===\n", fTd->fTid);
    int nworkers = propagator->fNthreads;
    for (int i = 0; i < nworkers; i++)
      wm->FeederQueue()->push(0);
    wm->TransportedQueue()->push(0);
    wm->Stop();
    //         sched_locker.StartOne(); // signal the scheduler who has to exit
    //         gbc_locker.StartOne();

    // WP
#ifdef USE_ROOT
    if (concurrentWrite) {
      file->Write();
    }
#endif

    wm->DoneQueue()->push(0);

    // Final reduction of counters
    propagator->fNsteps += fTd->fNsteps;
    propagator->fNsnext += fTd->fNsnext;
    propagator->fNphys += fTd->fNphys;
    propagator->fNmag += fTd->fNmag;
    propagator->fNsmall += fTd->fNsmall;
    propagator->fNcross += fTd->fNcross;

    if (wm->IsStopped()) wm->MergingServer()->Finish();

    if (concurrentWrite) {
      delete file;
    }

    return NULL;
  }else{
      // spawn feeder task
      while(runmgr->TryLock())
        ;
      //if(propagator->TryLock()){
      //  tbb::task::set_ref_count(2);
      //  TransportTask & transportTask = *new(tbb::task::allocate_child()) TransportTask( fTd );
      //  return & transportTask;
      //}
      tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
      FeederTask & feederTask = *new(cont.allocate_child()) FeederTask( fTd );
      return & feederTask;
  }
    return NULL;
}
