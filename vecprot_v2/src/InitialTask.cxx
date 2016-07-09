#include "InitialTask.h"
#include "ThreadData.h"
#include "FlowControllerTask.h"

#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif

using namespace Geant;

InitialTask::InitialTask () { }

InitialTask::~InitialTask () { }

tbb::task* InitialTask::execute ()
{

  GeantPropagator *propagator = GeantPropagator::Instance();
  WorkloadManager *wm = WorkloadManager::Instance();
  GeantScheduler *sch = wm->GetScheduler();
  ThreadData *threadData = ThreadData::Instance(propagator->fNthreads);

  tbb::task_scheduler_init init( propagator->fNthreads );
  int tid = wm->Instance()->ThreadId();
  Geant::Print("","=== Initial task %d created ===", tid);
  Geant::GeantTaskData *td = propagator->fThreadData[tid];
  td->fTid = tid;

  threadData->fPrioritizers[tid] = new GeantBasketMgr(sch, 0, 0, true);
  td->fBmgr = threadData->fPrioritizers[tid];
  threadData->fPrioritizers[tid]->SetThreshold(propagator->fNperBasket);
  threadData->fPrioritizers[tid]->SetFeederQueue(wm->FeederQueue());


  GeantFactoryStore* factoryStore = GeantFactoryStore::Instance();
  threadData->fMyhitFactories[tid] = factoryStore->GetFactory<MyHit>(16);

  bool concurrentWrite = GeantPropagator::Instance()->fConcurrentWrite && GeantPropagator::Instance()->fFillTree;
  if (concurrentWrite)
    {
      threadData->fFiles[tid] = new TThreadMergingFile("hits_output.root", wm->IOQueue(), "RECREATE");
      threadData->fTrees[tid] = new TTree("Tree","Simulation output");

      threadData->fTrees[tid]->Branch("hitblockoutput", "GeantBlock<MyHit>", &threadData->fData[tid]);

      // set factory to use thread-local queues
      threadData->fMyhitFactories[tid]->queue_per_thread = true;
    }
  tbb::task::set_ref_count(2);
  FlowControllerTask & flowControllerTask = *new(tbb::task::allocate_child()) FlowControllerTask();
  //tbb::task::spawn(flowControllerTask);
  return & flowControllerTask;
}
