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

  GeantBasketMgr *prioritizer = threadData->fPrioritizers[tid];
  prioritizer = new GeantBasketMgr(sch, 0, 0, true);
  td->fBmgr = prioritizer;
  prioritizer->SetThreshold(propagator->fNperBasket);
  prioritizer->SetFeederQueue(wm->FeederQueue());

  //TThread t;
  TThreadMergingFile* file = threadData->fFiles[tid];
  TTree *tree = threadData->fTrees[tid];
  GeantBlock<MyHit>* data = threadData->fData[tid];

  GeantFactoryStore* factoryStore = GeantFactoryStore::Instance();
  GeantFactory<MyHit> *myhitFactory = threadData->fMyhitFactories[tid];
  myhitFactory = factoryStore->GetFactory<MyHit>(16);

  bool concurrentWrite = GeantPropagator::Instance()->fConcurrentWrite && GeantPropagator::Instance()->fFillTree;
  if (concurrentWrite)
    {
      file = new TThreadMergingFile("hits_output.root", wm->IOQueue(), "RECREATE");
      tree = new TTree("Tree","Simulation output");

      tree->Branch("hitblockoutput", "GeantBlock<MyHit>", &data);

      // set factory to use thread-local queues
      myhitFactory->queue_per_thread = true;
    }

  FlowControllerTask & flowControllerTask = *new(tbb::task::allocate_child()) FlowControllerTask();
  tbb::task::spawn(flowControllerTask);
  return NULL;
}
