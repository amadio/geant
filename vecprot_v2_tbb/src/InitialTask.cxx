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

using namespace Geant;

tbb::task* InitialTask::execute ()
{
  GeantRunManager *runmgr = fPropagator->fRunMgr;
  WorkloadManager *wm = fPropagator->fWMgr;
  GeantScheduler *sch = wm->GetScheduler();
  ThreadData *threadData = ThreadData::Instance(runmgr->GetNthreadsTotal());

  tbb::task_scheduler_init init( runmgr->GetNthreadsTotal() );

  GeantTaskData *td = runmgr->GetTDManager()->GetTaskData();
  int tid = td->fTid;
  td->fPropagator = fPropagator;

  printf("=== Initial task %d (%d) created ===\n", tid, td->fTid);

  GeantBasketMgr *prioritizer = new GeantBasketMgr(fPropagator, sch, 0, 0, true);
  td->fBmgr = prioritizer;
  prioritizer->SetThreshold(fPropagator->fConfig->fNperBasket);
  prioritizer->SetFeederQueue(wm->FeederQueue());

  GeantFactoryStore* factoryStore = GeantFactoryStore::Instance();
  threadData->fMyhitFactories[tid] = factoryStore->GetFactory<MyHit>(16, runmgr->GetNthreadsTotal());

#ifdef USE_ROOT
  bool concurrentWrite = fPropagator->fConfig->fConcurrentWrite &&
                         fPropagator->fConfig->fFillTree;
  if (concurrentWrite)
    {
      threadData->fFiles[tid] = new TThreadMergingFile("hits_output.root", wm->IOQueue(), "RECREATE");
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
