#include "Geant/Error.h"

#include "FlowControllerTask.h"
#include "TransportTask.h"
#include "ThreadData.h"
#include "FeederTask.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "GeantEvent.h"
#include "GeantTaskData.h"
#include "GeantScheduler.h"
#include "GeantBasket.h"

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif
#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif
#include "TThread.h"
#include "GeantFactoryStore.h"
#include "TThreadMergingFile.h"

#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif


FlowControllerTask::FlowControllerTask (bool starting): fStarting(starting) { }

FlowControllerTask::~FlowControllerTask () { }

tbb::task* FlowControllerTask::execute ()
{
  WorkloadManager *wm = WorkloadManager::Instance();
  GeantPropagator *propagator = GeantPropagator::Instance();
  int tid = wm->ThreadId();

  // verification, just in case
  if(tid>=propagator->fNthreads)
    return NULL;

  Geant::GeantTaskData *td = propagator->fThreadData[tid];

  if(fStarting){
    while(propagator->TryLock())
      ;

    int didFeeder;
    tbb::task::set_ref_count(2);
    FeederTask & feederTask = *new(tbb::task::allocate_child()) FeederTask(&didFeeder);
    return & feederTask;
  }

  GeantScheduler *sch = wm->GetScheduler();
  ThreadData *threadData = ThreadData::Instance(propagator->fNthreads);
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();

  GeantBasketMgr *prioritizer = threadData->fPrioritizers[tid];
  Geant::TThreadMergingFile* file = threadData->fFiles[tid];
  TTree *tree = threadData->fTrees[tid];
  GeantBlock<MyHit>* data = threadData->fData[tid];
  int nbaskets = 0;
  int nworkers = propagator->fNthreads;
  int *waiting = wm->GetWaiting();

  bool concurrentWrite = propagator->fConcurrentWrite && propagator->fFillTree;

  if (propagator->TransportCompleted()) {
    printf("=== Exit thread %d from Flow Controller  ===\n", tid);
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

    delete prioritizer;
    // Final reduction of counters
    propagator->fNsteps += td->fNsteps;
    propagator->fNsnext += td->fNsnext;
    propagator->fNphys += td->fNphys;
    propagator->fNmag += td->fNmag;
    propagator->fNsmall += td->fNsmall;
    propagator->fNcross += td->fNcross;

    if (wm->IsStopped()) wm->MergingServer()->Finish();

    if (concurrentWrite) {
      delete file;
    }
    return NULL;
  }else{

      // spawn feeder task
      while(propagator->TryLock())
        ;

      int didFeeder;
      tbb::task::set_ref_count(2);
      FeederTask & feederTask = *new(tbb::task::allocate_child()) FeederTask(&didFeeder);
      return & feederTask;
  }
    return NULL;
}
