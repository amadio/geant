#include "TTree.h"
#include "GeantPropagator.h"
#include "WorkloadManager.h"
#include "GeantScheduler.h"
#include "GeantTaskData.h"
#include "GeantBasket.h"
#include "MyHit.h"
#include "TThreadMergingFile.h"

#include "InitialTask.h"

#include "tbb/task_scheduler_init.h"

tbb::task* InitialTask::execute ()
{
  GeantPropagator *propagator = GeantPropagator::Instance();
  WorkloadManager *wm = WorkloadManager::Instance();
  GeantScheduler *sch = wm->GetScheduler();
  ThreadData *threadData = ThreadData::Instance(propagator->fNthreads);

  tbb::task_scheduler_init init( propagator->fNthreads );

  int tid = wm->Instance()->ThreadId(); // only for checking thread id with initial task
  Geant::GeantTaskData *td = propagator->fThreadData[tid];

  printf("=== Initial task %d (%d) created ===\n", tid, td->fTid);

  threadData->fPrioritizers[tid] = new GeantBasketMgr(sch, 0, 0, true);
  td->fBmgr = threadData->fPrioritizers[tid];
  threadData->fPrioritizers[tid]->SetThreshold(propagator->fNperBasket);
  threadData->fPrioritizers[tid]->SetFeederQueue(wm->FeederQueue());

  GeantFactoryStore* factoryStore = GeantFactoryStore::Instance();
  threadData->fMyhitFactories[tid] = factoryStore->GetFactory<MyHit>(16);

  bool concurrentWrite = GeantPropagator::Instance()->fConcurrentWrite && GeantPropagator::Instance()->fFillTree;
  if (concurrentWrite)
    {
      threadData->fFiles[tid] = new Geant::TThreadMergingFile("hits_output.root", wm->IOQueue(), "RECREATE");
      threadData->fTrees[tid] = new TTree("Tree","Simulation output");

      threadData->fTrees[tid]->Branch("hitblockoutput", "GeantBlock<MyHit>", &threadData->fData[tid]);

      // set factory to use thread-local queues
      threadData->fMyhitFactories[tid]->queue_per_thread = true;
    }

  tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
  FlowControllerTask & flowControllerTask = *new(cont.allocate_child()) FlowControllerTask(td, true);
  return & flowControllerTask;
}
