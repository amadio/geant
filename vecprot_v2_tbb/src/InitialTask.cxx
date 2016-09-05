#include "InitialTask.h"

#include "GeantRunManager.h"
#include "GeantPropagator.h"
#include "WorkloadManager.h"
#include "GeantScheduler.h"
#include "GeantTaskData.h"
#include "GeantBasket.h"
#include "GeantFactory.h"
#include "MyHit.h"

#ifdef USE_ROOT
#include "TTree.h"
#include "TThreadMergingFile.h"
#endif

#include "tbb/task_scheduler_init.h"
#include "ThreadData.h"
#include "FlowControllerTask.h"

tbb::task* InitialTask::execute ()
{
  GeantRunManager *runmgr = fPropagator->fRunMgr;
  WorkloadManager *wm = fPropagator->fWMgr;
  GeantScheduler *sch = wm->GetScheduler();
  ThreadData *threadData = ThreadData::Instance(runmgr->GetNthreadsTotal());

  tbb::task_scheduler_init init( runmgr->GetNthreadsTotal() );

  int tid = runmgr->GetTaskId();
  Geant::GeantTaskData *td = runmgr->GetTaskData(tid);

  printf("=== Initial task %d (%d) created ===\n", tid, td->fTid);

  threadData->fPrioritizers[tid] = new GeantBasketMgr(fPropagator, sch, 0, 0, true);
  td->fBmgr = threadData->fPrioritizers[tid];
  threadData->fPrioritizers[tid]->SetThreshold(fPropagator->fConfig->fNperBasket);
  threadData->fPrioritizers[tid]->SetFeederQueue(wm->FeederQueue());

  GeantFactoryStore* factoryStore = GeantFactoryStore::Instance();
  threadData->fMyhitFactories[tid] = factoryStore->GetFactory<MyHit>(16, runmgr->GetNthreadsTotal());

#ifdef USE_ROOT
  bool concurrentWrite = fPropagator->fConfig->fConcurrentWrite &&
                         fPropagator->fConfig->fFillTree;
  if (concurrentWrite)
    {
      threadData->fFiles[tid] = new Geant::TThreadMergingFile("hits_output.root", wm->IOQueue(), "RECREATE");
      threadData->fTrees[tid] = new TTree("Tree","Simulation output");

      threadData->fTrees[tid]->Branch("hitblockoutput", "GeantBlock<MyHit>", &threadData->fData[tid]);

      // set factory to use thread-local queues
      threadData->fMyhitFactories[tid]->queue_per_thread = true;
    }
#endif

  tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
  FlowControllerTask & flowControllerTask = *new(cont.allocate_child()) FlowControllerTask(td, true);
  return & flowControllerTask;
}
