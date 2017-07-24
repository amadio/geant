#include "FlowControllerTask.h"

#include "ThreadData.h"
#include "FeederTask.h"
#include "TransportTask.h"
#include "GeantRunManager.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"

#ifdef USE_ROOT
#include "TThreadMergingFile.h"
#endif

#include "tbb/task_scheduler_init.h"

using namespace Geant;

FlowControllerTask::FlowControllerTask (Geant::GeantTaskData *td, bool starting, bool forcedStop)
  :fTd(td), fStarting(starting), fForcedStop(forcedStop) { }

FlowControllerTask::~FlowControllerTask () { }

tbb::task* FlowControllerTask::SpawnFeederTask()
{
  tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
  FeederTask & feederTask = *new(cont.allocate_child()) FeederTask( fTd, fStarting );
  return & feederTask;
}

tbb::task* FlowControllerTask::SpawnTransportTask()
{
  tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
  TransportTask & transportTask = *new(cont.allocate_child()) TransportTask( fTd, fStarting );
  return & transportTask;
}

tbb::task* FlowControllerTask::execute ()
{
  GeantPropagator *propagator = fTd->fPropagator;
  GeantRunManager *runmgr = propagator->fRunMgr;
  WorkloadManager *wm = propagator->fWMgr;
  bool multiPropagator = runmgr->GetNpropagators() > 1;

  //printf("=== %d Flow Controller  ===\n", fTd->fTid);
  if(fStarting) {
//    if (runmgr->GetFedPropagator() == propagator)
      return SpawnTransportTask(); 

//    return SpawnFeederTask();
  }

#ifdef USE_ROOT
  ThreadData *threadData = ThreadData::Instance(runmgr->GetNthreadsTotal());
  Geant::TThreadMergingFile* file = threadData->fFiles[fTd->fTid];
  bool concurrentWrite = propagator->fConfig->fConcurrentWrite &&
                         propagator->fConfig->fFillTree;
#endif

  if (fForcedStop || runmgr->TransportCompleted()) {
    //finish tasks
    //tbb::task::destroy(*tbb::task::parent());
    printf("=== Exit thread %d from Flow Controller  ===\n", fTd->fTid);
    // Stop other idle propagators
    if (multiPropagator) {
      GeantPropagator *idle = propagator->fRunMgr->GetIdlePropagator();
      if (idle) idle->StopTransport();
    }
    int nworkers = propagator->fNthreads;
    if (!fForcedStop) {
      // Stop this propagator
      for (int i = 0; i < nworkers; i++)
        wm->FeederQueue()->push(0);
    }
    wm->TransportedQueue()->push(0);
    wm->Stop();
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

#ifdef USE_ROOT
    if (wm->IsStopped()) wm->MergingServer()->Finish();

    if (concurrentWrite) {
      delete file;
    }
#endif
  } else {
    // spawn feeder task
//    if (runmgr->IsFeeding(propagator))
      return SpawnTransportTask(); 
//    while(runmgr->IsFeeding(propagator))
//      ;
    //if(propagator->TryLock()){
    //  tbb::task::set_ref_count(2);
    //  TransportTask & transportTask = *new(tbb::task::allocate_child()) TransportTask( fTd );
    //  return & transportTask;
    //}
//    return SpawnFeederTask();
  }
  return NULL;
}
