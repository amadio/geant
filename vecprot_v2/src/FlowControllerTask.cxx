#include "Geant/Error.h"

#include "FlowControllerTask.h"
#include "ThreadData.h"
#include "FeederTask.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "GeantEvent.h"
#include "GeantTaskData.h"

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

using namespace Geant;

FlowControllerTask::FlowControllerTask () { }

FlowControllerTask::~FlowControllerTask () { }

tbb::task* FlowControllerTask::execute ()
{
  GeantPropagator *propagator = GeantPropagator::Instance();
  WorkloadManager *wm = WorkloadManager::Instance();
  GeantScheduler *sch = wm->GetScheduler();
    ThreadData *threadData = ThreadData::Instance(propagator->fNthreads);

  int tid = wm->Instance()->ThreadId();
  GeantTaskData *td = propagator->fThreadData[tid];
  GeantBasketMgr *prioritizer = threadData->fPrioritizers[tid];
  TThreadMergingFile* file = threadData->fFiles[tid];
  TTree *tree = threadData->fTrees[tid];
  GeantBlock<MyHit>* data = threadData->fData[tid];


  Geant::Print("","=== FlowControllerTask task %d created ===", tid);

  //int tid = wm->Instance()->ThreadId();
  //Geant::GeantTaskData *td = propagator->fThreadData[tid];

  bool concurrentWrite = propagator->fConcurrentWrite && propagator->fFillTree;

  // Fire garbage collection if starving
  /*
  if ((*nbaskets < 1) && (!propagator->IsFeeding())) {
    printf("garbage collector \n");
    sch->GarbageCollect(fTd);
    *fNgcoll++;

    // Too many garbage collections - enter priority mode
    if ((*fNgcoll > 5) && (wm->GetNworking() <= 1)) {
      fNgcoll = 0;
      for (int slot = 0; slot < propagator->fNevents; slot++)
        if (propagator->fEvents[slot]->Prioritize())
          propagator->fPriorityEvents++;
      while ((!sch->GarbageCollect(fTd, true)) &&
             (wm->FeederQueue()->size_async() == 0) &&
             (!*nbaskets) &&
             (!(fPrioritizer->HasTracks())))
        ;
    }
  }else{*/
    // Check exit condition
    if (propagator->TransportCompleted()) {
      printf("exit...\n");
      int nworkers = propagator->fNthreads;
      for (int i = 0; i < nworkers; i++)
        wm->FeederQueue()->push(0);
      wm->TransportedQueue()->push(0);
      wm->Stop();
      //         sched_locker.StartOne(); // signal the scheduler who has to exit
      //         gbc_locker.StartOne();


      // Here spawn I/O task
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
    }else{
      Geant::Print("","=== creating feeder task %d  ===", tid);
      // spawn feederQ
      //bool didFeeder = propagator.Feeder(&td);
      int didFeeder;
      FeederTask & feederTask = *new(tbb::task::allocate_child()) FeederTask(&didFeeder);
      tbb::task::spawn(feederTask);
    }

  //}
    return NULL;
}
