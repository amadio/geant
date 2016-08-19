#include "FlowControllerTask.h"

#include "ThreadData.h"
#include "FeederTask.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "TThreadMergingFile.h"

#include "tbb/task_scheduler_init.h"


FlowControllerTask::FlowControllerTask (Geant::GeantTaskData *td, bool starting): fTd(td), fStarting(starting) { }

FlowControllerTask::~FlowControllerTask () { }

tbb::task* FlowControllerTask::execute ()
{
  WorkloadManager *wm = WorkloadManager::Instance();
  GeantPropagator *propagator = GeantPropagator::Instance();

  //printf("=== %d Flow Controller  ===\n", fTd->fTid);
  if(fStarting){
    while(propagator->TryLock())
      ;

    tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
    FeederTask & feederTask = *new(cont.allocate_child()) FeederTask( fTd );
    return & feederTask;
  }

  ThreadData *threadData = ThreadData::Instance(propagator->fNthreads);
  Geant::TThreadMergingFile* file = threadData->fFiles[fTd->fTid];
  bool concurrentWrite = propagator->fConcurrentWrite && propagator->fFillTree;

  if (propagator->TransportCompleted()) {
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
    if (concurrentWrite) {
      file->Write();
    }

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
      while(propagator->TryLock())
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
