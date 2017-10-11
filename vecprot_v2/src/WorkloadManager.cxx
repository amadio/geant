#include "WorkloadManager.h"

#include "Geant/Error.h"

#ifdef USE_ROOT
#include "TApplication.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TBits.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#endif

#include "GeantRunManager.h"
#include "GeantEventServer.h"
#include "GeantTrackVec.h"
#include "GeantBasket.h"
#include "GeantOutput.h"
#include "GeantTaskData.h"
#include "PhysicsInterface.h"
#include "PhysicsProcessOld.h"
#include "GeantScheduler.h"
#include "GeantEvent.h"
#include "GeantVApplication.h"
#include "GeantVTaskMgr.h"
#include "StackLikeBuffer.h"
#include "SimulationStage.h"
#include "TrackStat.h"
#include "LocalityManager.h"
#include "TrackManager.h"
#if USE_VECGEOM_NAVIGATOR
#include "base/TLS.h"
#include "management/GeoManager.h"
#include "base/Stopwatch.h"
#else
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#endif
#include "TaskBroker.h"

// added by WP for output handling
#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif
#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif
#ifdef USE_ROOT
#include "TThreadMergingFile.h"
#endif
#include "GeantFactoryStore.h"

using std::max;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
WorkloadManager::WorkloadManager(int nthreads, GeantPropagator* prop)
    : fPropagator(prop), fNthreads(nthreads), fNbaskets(0), fBasketGeneration(0), fNbasketgen(0), fNidle(nthreads),
      fNqueued(0), fBtogo(0), fSchId(nthreads), fStarted(false), fStopped(false), fFeederQ(0), fTransportedQ(0),
      fDoneQ(0), fListThreads(), fFlushed(false), fFilling(false), fScheduler(0), fBroker(0), fSchLocker(), fGbcLocker()
#ifdef USE_ROOT
    , fOutputIO(0)
#endif
  {
  // Private constructor.
  fFeederQ = new Geant::priority_queue<GeantBasket *>(1 << 16, nthreads);
  fTransportedQ = new Geant::priority_queue<GeantBasket *>(1 << 16, nthreads);
  fDoneQ = new Geant::priority_queue<GeantBasket *>(1 << 10, nthreads);
  fScheduler = new GeantScheduler();
  ReleaseShareLock();

#ifdef USE_ROOT
  fOutputIO = new dcqueue<TBufferFile*>;
  fMergingServer = new TThreadMergingServer(fOutputIO);
#endif
}

//______________________________________________________________________________
WorkloadManager::~WorkloadManager() {
  // Destructor.
  delete fFeederQ;
  delete fTransportedQ;
  delete fDoneQ;
  delete fScheduler;
#ifdef USE_ROOT
  delete fOutputIO;
#endif
}

//______________________________________________________________________________
int WorkloadManager::ThreadId() {
#if USE_VECGEOM_NAVIGATOR
  return BaseTLS::ThreadId();
#else
  return TGeoManager::ThreadId();
#endif
}

//______________________________________________________________________________
void WorkloadManager::CreateBaskets(GeantPropagator* prop) {
  // Create the array of baskets
  VolumePath_t *blueprint = 0;
  int maxdepth = prop->fConfig->fMaxDepth;
  blueprint = VolumePath_t::MakeInstance(maxdepth);
  fScheduler->CreateBaskets(prop);
  VolumePath_t::ReleaseInstance(blueprint);

  if (fBroker)
    fBroker->CreateBaskets(prop);
}

//______________________________________________________________________________
WorkloadManager *WorkloadManager::NewInstance(GeantPropagator *prop, int nthreads) {
  // Return singleton instance.
  if (!nthreads || !prop) {
    Geant::Error("WorkloadManager::NewInstance", "%s", "You should provide number of threads.");
    return 0;
  }
  return new WorkloadManager(nthreads,prop);
}

//______________________________________________________________________________
void WorkloadManager::Print(const char *) const {
  //
}

void WorkloadManager::SetTaskBroker(TaskBroker *broker) {
  // Register the broker if it is valid.
  if (broker && broker->IsValid())
    fBroker = broker;
}

#ifdef USE_VECGEOM_NAVIGATOR
//______________________________________________________________________________
bool WorkloadManager::LoadGeometry(vecgeom::VPlacedVolume const *const volume) {
  /**
   * @brief Tell the task broker(s) to load the geometry.
   *
   * @param Volume to load
   */
  if (fBroker)
    return fBroker->UploadGeometry(volume);
  return true;
}
#endif

//______________________________________________________________________________
bool WorkloadManager::StartTasks(GeantVTaskMgr *taskmgr) {
  GeantPropagator *prop=fPropagator;
  // Start the threads
  fStarted = true;
  if (!fListThreads.empty())
    return false;
  int ith = 0;
  if (fBroker) {
     if (fBroker->GetNstream() > (unsigned int)fNthreads) {
       Geant::Fatal("StartThreads", "The task broker is using too many threads (%d out of %d)", fBroker->GetNstream(),
                    fNthreads);
       return false;
    }
    Geant::Info("StartThreads", "Running with a coprocessor broker (using %d threads).",fBroker->GetNstream()+1);
    fListThreads.emplace_back(WorkloadManager::TransportTracksCoprocessor, prop, fBroker);
    ith += fBroker->GetNstream() + 1;
    if (ith == fNthreads && fBroker->IsSelective()) {
       Geant::Fatal("WorkloadManager::StartThreads","All %d threads are used by the coprocessor broker but it can only process a subset of particles.",fNthreads);
       return false;
    }
  }

  // Start CPU transport threads (static mode)
  if (!taskmgr) {
    for (; ith < fNthreads; ith++) {
      if ( fPropagator->fConfig->fUseV3 )
        fListThreads.emplace_back(WorkloadManager::TransportTracksV3, prop);
      else
        fListThreads.emplace_back(WorkloadManager::TransportTracks, prop);
    }
  }

  // *** THIS IS TO BE MOVED IN THE RUN MANAGER ***

  // Start output thread
  if (prop->fConfig->fFillTree) {
    fListThreads.emplace_back(WorkloadManager::OutputThread, prop);
  }
  // Start monitoring thread
  if (prop->fConfig->fUseMonitoring) {
    //fListThreads.emplace_back(WorkloadManager::StartROOTApplication);
    fListThreads.emplace_back(WorkloadManager::MonitoringThread, prop);
  }
  // Start garbage collector
  if (prop->fConfig->fMaxRes > 0) {
    fListThreads.emplace_back(WorkloadManager::GarbageCollectorThread, prop);
  }
  if (taskmgr) {
    Printf("%s","=== TBB Task Mode ====");
    return taskmgr->Initialize(fNthreads, prop);
  } else {
    Printf("%s","=== Thread Mode ====");
  }
  return true;
}

//______________________________________________________________________________
void WorkloadManager::JoinThreads() {
  //
  StopTransportThreads();
  for (auto &t : fListThreads) {
    t.join();
  }
}

//______________________________________________________________________________
void WorkloadManager::StopTransportThreads() {
  int tojoin = fNthreads;
  for (int ith = 0; ith < tojoin; ith++)
    fFeederQ->push_force(0);
  fStopped = true;
}

//______________________________________________________________________________
void WorkloadManager::WaitWorkers() {
  // Waiting point for the main thread until work gets done.
  fFilling = false;
  int ntowait = fNthreads;
  if (fBroker)
    ntowait -= fBroker->GetNstream();
  GeantBasket *signal;
  while (ntowait) {
    fDoneQ->wait_and_pop(signal);
    ntowait--;
    Geant::Print("", "=== %d workers finished", fNthreads - ntowait);
  }
  //   fBasketGeneration++;
}

//______________________________________________________________________________
void *WorkloadManager::MainScheduler(void *) {
  // Garbage collector thread, called by a single thread.
  Geant::Print("","%s","=== Scheduler: stopping threads and exiting ===");
  return 0;
}

//______________________________________________________________________________
static inline void MaybeCleanupBaskets(GeantTaskData *td) {
  if (td->NeedsToClean())
    td->CleanBaskets(0);
  else {
    // Check if there are at least 10 free baskets for this thread
    if (td->fPool.size() < 10) {
       td->fBmgr->CreateEmptyBaskets(10, td);
    }
  }
}

//______________________________________________________________________________
static inline void MaybeCleanupBaskets(GeantTaskData *td, GeantBasket *basket) {
  if (td->NeedsToClean())
    td->CleanBaskets(0);
  else {
    // Check if there are at least 10 free baskets for this thread
    if (td->fPool.size() < 10) {
      if (basket->GetBasketMgr())
        basket->GetBasketMgr()->CreateEmptyBaskets(2, td);
      else
        td->fBmgr->CreateEmptyBaskets(10, td);
    }
  }
}

//______________________________________________________________________________
/** @brief Call Feeder (if needed) and check exit condition. */
WorkloadManager::FeederResult WorkloadManager::CheckFeederAndExit() {

//  if (!prioritizer.HasTracks() && (propagator.fRunMgr->GetNpriority() || propagator.GetNworking() == 1)) {
//    bool didFeeder = propagator.fRunMgr->Feeder(&td) > 0;

  // Check exit condition
  if (fPropagator->fRunMgr->TransportCompleted()) {
    int nworkers = fPropagator->fNthreads;
    for (int i = 0; i < nworkers; i++)
      FeederQueue()->push_force(0);
    TransportedQueue()->push_force(0);
    Stop();
    return FeederResult::kStop;
  }
  return FeederResult::kNone;
}

//______________________________________________________________________________
int WorkloadManager::ShareBaskets(WorkloadManager *other)
{
// Share some baskets with a different workload manager
  // Check transport queue status and share only if enough baskets
  if (TryShareLock())
    return 0;
  unsigned int status = fFeederQ->status();
  if (status < 2) {
    ReleaseShareLock();
    return 0;
  }
  GeantBasket *basket = nullptr;
  int ntoshare = Math::Min<int>(fFeederQ->size() - fNthreads, fNthreads);
  ntoshare = Math::Max<int>(ntoshare, 0);
  int nshared = 0;
  for (auto i=0; i<ntoshare; ++i) {
    if (fFeederQ->try_pop(basket)) {
      other->FeederQueue()->push_force(basket);
      nshared++;
    }
  }
  ReleaseShareLock();
  return nshared;
}
//______________________________________________________________________________
void WorkloadManager::TransportTracksV3(GeantPropagator *prop) {

//  int nstacked = 0;
//  int nprioritized = 0;
//  int ninjected = 0;
  bool useNuma = prop->fConfig->fUseNuma;
  // Enforce locality by pinning the thread to the next core according to the chosen policy.
  int node = -1;
  LocalityManager *loc_mgr = LocalityManager::Instance();
  if (useNuma) node = loc_mgr->GetPolicy().AllocateNextThread(prop->fNuma);
  int cpu = useNuma ? NumaUtils::GetCpuBinding() : -1;
//  if (node < 0) node = 0;
  GeantPropagator *propagator = prop;
  GeantRunManager *runmgr = prop->fRunMgr;
  Geant::GeantTaskData *td = runmgr->GetTDManager()->GetTaskData();
  td->AttachPropagator(prop, node);
  int tid = td->fTid;

  if (useNuma) {
    Geant::Print("","=== Worker thread %d created for propagator %p on NUMA node %d CPU %d ===",
                 tid, prop, node, cpu);    
    int membind = NumaUtils::NumaNodeAddr(td->fShuttleBasket->Tracks().data());
    if (node != membind)
      Geant::Print("","=== Thread #d: Wrong memory binding");
  } else {
    Geant::Print("","=== Worker thread %d created for propagator %p ===", tid, prop);
  }

  GeantEventServer *evserv = runmgr->GetEventServer();

  // IO handling
  #ifdef USE_ROOT
  bool concurrentWrite = prop->fConfig->fConcurrentWrite && prop->fConfig->fFillTree;
//  int treeSizeWriteThreshold = td->fPropagator->fConfig->fTreeSizeWriteThreshold;

  GeantFactoryStore* factoryStore = GeantFactoryStore::Instance();
  GeantFactory<MyHit> *myhitFactory = factoryStore->GetFactory<MyHit>(16,runmgr->GetNthreadsTotal());

  TThreadMergingFile* file=0;
  TTree *tree=0;
  GeantBlock<MyHit>* data=0;

  if (concurrentWrite)
    {
      file = new TThreadMergingFile("hits_output.root", propagator->fWMgr->IOQueue(), "RECREATE");
      tree = new TTree("Tree","Simulation output");

      tree->Branch("hitblockoutput", "GeantBlock<MyHit>", &data);

      // set factory to use thread-local queues
      myhitFactory->queue_per_thread = true;
    }
  #endif

#ifndef USE_VECGEOM_NAVIGATOR
  // If we use ROOT make sure we have a navigator here
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif

/*** STEPPING LOOP ***/
  bool flush = false;
  while (1) {
    // Check if simulation has finished
    auto feedres = PreloadTracksForStep(td); 
    if (feedres == FeederResult::kStop) break;
    if (flush) {
      if ((feedres == FeederResult::kNone) | (feedres == FeederResult::kError)) {
        if (!evserv->EventsServed())
          Geant::Warning("","=== Task %d exited due to missing workload", tid);
        break;
      }
    }
    flush = (feedres == FeederResult::kNone) | (feedres == FeederResult::kError);
    SteppingLoop(td, flush);
  }
  // WP
  #ifdef USE_ROOT
  if(concurrentWrite)
    {
      file->Write();
    }
  #endif
  prop->fWMgr->DoneQueue()->push_force(nullptr);
  // Final reduction of counters
  propagator->fNsteps += td->fNsteps;
  propagator->fNsnext += td->fNsnext;
  propagator->fNphys += td->fNphys;
  propagator->fNmag += td->fNmag;
  propagator->fNsmall += td->fNsmall;
  propagator->fNcross += td->fNcross;
  propagator->fNpushed += td->fNpushed;
  propagator->fNkilled += td->fNkilled;

  // If transport is completed, make send the signal to the run manager
  if (runmgr->TransportCompleted())
    runmgr->StopTransport();
  if (useNuma) {
    int cpuexit = NumaUtils::GetCpuBinding();
    if (cpuexit != cpu)
      Geant::Print("","=== OS migrated worker %d from cpu #%d to cpu#%d", tid, cpu, cpuexit);
  }
  runmgr->GetTDManager()->ReleaseTaskData(td);
  Geant::Print("","=== Thread %d: exiting ===", tid);

  #ifdef USE_ROOT
  if (prop->fWMgr->IsStopped()) prop->fWMgr->MergingServer()->Finish();

  if (concurrentWrite) {
    delete file;
  }
  #endif
}

//______________________________________________________________________________
WorkloadManager::FeederResult WorkloadManager::PreloadTracksForStep(GeantTaskData *td) {
 // The method will apply a policy to inject tracks into the first stepping
 // stage buffer. The policy has to be exhaustively described...
  auto feedres = td->fPropagator->fWMgr->CheckFeederAndExit();
  if (feedres == FeederResult::kStop)
     return feedres;

  // Take tracks from the event server
  int ninjected = 0;
  unsigned int error = 0;
  GeantEventServer *evserv = td->fPropagator->fRunMgr->GetEventServer();
  if (evserv->HasTracks()) {
    // In the initial phase we distribute a fair share of baskets to all propagators
    if (!evserv->IsInitialPhase() ||
        td->fPropagator->fNbfeed < td->fPropagator->fRunMgr->GetInitialShare())
      ninjected = evserv->FillStackBuffer(td->fStackBuffer, td->fPropagator->fConfig->fNperBasket, error);
  }
  // td->fStat->AddTracks(ninjected);
  if (ninjected)
    return FeederResult::kWork;

  if (error)
    return FeederResult::kError;

  return FeederResult::kNone;
}

//______________________________________________________________________________
int WorkloadManager::FlushOneLane(GeantTaskData *td)
{
// Flush a single track lane from the stack-like buffer into the first stage.
  // Check the stack buffer and flush priority events first
  int ninjected = 0;
  int maxspill = td->fPropagator->fConfig->fNmaxBuffSpill;
  if ( td->fStackBuffer->IsPrioritized())
    ninjected = td->fStackBuffer->FlushPriorityLane();
  if (ninjected) return ninjected;

  // Inject the last lane in the buffer
  int nlane = td->fStackBuffer->FlushLastLane();
  ninjected += nlane;
  while (nlane && (ninjected < maxspill)) {
    nlane = td->fStackBuffer->FlushLastLane();
    ninjected += nlane;
  }
  return ( ninjected );
}

//______________________________________________________________________________
int WorkloadManager::SteppingLoop(GeantTaskData *td, bool flush)
{
// The main stepping loop over simulation stages.
//  static int count = 0;
  constexpr int nstages = int(ESimulationStage::kSteppingActionsStage) + 1;
  bool flushed = !flush;
  int nprocessed = 0;
  int ninput = 0;
  int istage = 0;
  while ( FlushOneLane(td) || !flushed ) {
    while (1) {
//      count++;
      int nstart = td->fStageBuffers[istage]->size();
/*
      for (auto track : td->fStageBuffers[istage]->Tracks()) {
        if (track->fEvent == 0 && track->fParticle == 0) {
          td->InspectStages(istage);
          track->Print("");
          break;
        }
      }
*/
      ninput += nstart;
      if ( nstart || !flushed ) {
        if (flush)
          nprocessed += td->fPropagator->fStages[istage]->FlushAndProcess(td);
        else
          nprocessed += td->fPropagator->fStages[istage]->Process(td);
      }
      istage = (istage + 1) % nstages;
      if (istage == 0) {
        if (flush) flushed = true;
        if (ninput == 0) break;
        ninput = 0;
      }
//      assert(td->fStat->CountBalance() == 0);
    }
  }
  return nprocessed; // useless instruction intended for dummy compilers
}

//______________________________________________________________________________
void *WorkloadManager::TransportTracks(GeantPropagator *prop) {
  // Thread propagating all tracks from a basket.
  //      char slist[256];
  //      TString sslist;
  //   const int max_idle = 1;
  //   int indmin, indmax;
  static std::atomic<int> counter(0);
  int ntotnext, nbaskets;
  int ntotransport;
  int ncross = 0;
  int nextra_at_rest = 0;
  int generation = 0;
//  int ninjected = 0;
  int nnew = 0;
  int ntot = 0;
  int nkilled = 0;
  int nphys = 0;
  int nout = 0;
  int ngcoll = 0;
  unsigned int error = 0;
  GeantBasket *basket = 0;
  GeantPropagator *propagator = prop;
  int basket_size = propagator->fConfig->fNperBasket;
  // Not NUMA aware version
  propagator->SetNuma(0);
  GeantRunManager *runmgr = prop->fRunMgr;
  Geant::GeantTaskData *td = runmgr->GetTDManager()->GetTaskData();
  int tid = td->fTid;
  td->AttachPropagator(prop, 0);
  Geant::Print("","=== Worker thread %d created for propagator %p ===", tid, prop);
  td->fPropagator = prop;
//  td->fStackBuffer = new StackLikeBuffer(propagator->fConfig->fNstackLanes, td);
//  td->fStackBuffer->SetStageBuffer(td->fStageBuffers[0]);
  int nworkers = propagator->fNthreads;
  WorkloadManager *wm = propagator->fWMgr;
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  GeantScheduler *sch = wm->GetScheduler();
  int *nvect = sch->GetNvect();
  GeantBasketMgr *prioritizer = new GeantBasketMgr(prop,sch, 0, 0, true);
  td->fBmgr = prioritizer;
  prioritizer->SetThreshold(basket_size);
  prioritizer->SetFeederQueue(feederQ);

  GeantEventServer *evserv = runmgr->GetEventServer();
  int bindex = evserv->GetBindex();
  GeantBasket *bserv = sch->GetBasketManagers()[bindex]->GetNextBasket(td);
  bserv->SetThreshold(basket_size);
  td->fImported = bserv;

  bool firstTime = true;
  bool multiPropagator = runmgr->GetNpropagators() > 1;
  GeantPropagator *idle = nullptr;

  // IO handling
  #ifdef USE_ROOT
  bool concurrentWrite = td->fPropagator->fConfig->fConcurrentWrite && td->fPropagator->fConfig->fFillTree;
  int treeSizeWriteThreshold = td->fPropagator->fConfig->fTreeSizeWriteThreshold;

  GeantFactoryStore* factoryStore = GeantFactoryStore::Instance();
  GeantFactory<MyHit> *myhitFactory = factoryStore->GetFactory<MyHit>(16,runmgr->GetNthreadsTotal());

  TThreadMergingFile* file=0;
  TTree *tree=0;
  GeantBlock<MyHit>* data=0;

  if (concurrentWrite)
    {
      file = new TThreadMergingFile("hits_output.root", wm->IOQueue(), "RECREATE");
      tree = new TTree("Tree","Simulation output");

      tree->Branch("hitblockoutput", "GeantBlock<MyHit>", &data);

      // set factory to use thread-local queues
      myhitFactory->queue_per_thread = true;
    }
  #endif

  // Start the feeder for this propagator
/*
  while (runmgr->GetFedPropagator() != propagator) {
    int nb0 = runmgr->Feeder(td);
    if (nb0 > 0) ninjected += nb0;
  }
  // The first injection must produce some baskets
  if (!ninjected) sch->GarbageCollect(td, true);
*/

  // Activate events in the server
  //evserv->ActivateEvents();

  Material_t *mat = 0;
#ifndef USE_VECGEOM_NAVIGATOR
  // If we use ROOT make sure we have a navigator here
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif
  //   int iev[500], itrack[500];
  // TGeoBranchArray *crt[500], *nxt[500];
  while (1) {
    // Call the feeder if in priority mode
    auto feedres = wm->CheckFeederAndExit();
    if (feedres == FeederResult::kStop) {
       break;
    }
    // Collect info about the queue
    nbaskets = feederQ->size_async();
    if (nbaskets > nworkers)
      ngcoll = 0;

    // Fire garbage collection if starving
//    if (!firstTime && (nbaskets < 1) && (!runmgr->IsFeeding(propagator))) {
    if (!firstTime && (nbaskets < 1) && (!evserv->HasTracks())) {
      sch->GarbageCollect(td);
     ngcoll++;
    }
    // Too many garbage collections - enter priority mode
    if ((ngcoll > 5) && (propagator->GetNworking() <= 1)) {
      ngcoll = 0;
      // std::atomic_int &priority_events = propagator->fRunMgr->GetPriorityEvents();
      // Max number of prioritized events should be configurable

      while (!sch->GarbageCollect(td, true) &&
             !feederQ->size_async() &&
             !basket &&
             !prioritizer->HasTracks() &&
             !evserv->HasTracks())
        ;
    }
    // Check if the current basket is reused or we need a new one
    if (!basket) {
      // If prioritizer has work, just do it
      if (prioritizer->HasTracks()) {
        basket = prioritizer->GetBasketForTransport(td);
        ngcoll = 0;
      } else {
        // Get basket from the generator
        if (evserv->HasTracks()) {
          // In the initial phase we distribute a fair share of baskets to all propagators
          if (!evserv->IsInitialPhase() ||
              propagator->fNbfeed < runmgr->GetInitialShare()) {
            ntotransport = evserv->FillBasket(bserv->GetInputTracks(), basket_size, error);
            if (ntotransport) basket = bserv;
            else sch->GarbageCollect(td);
          }
        }
        // Try to get from work queue without waiting
        if (!basket)
          wm->FeederQueue()->try_pop(basket);
        if (!basket) {
          // Take next basket from queue
          if (propagator->fCompleted || wm->fStopped)
            break;
          propagator->fNidle++;
          wm->FeederQueue()->wait_and_pop(basket);
          propagator->fNidle--;
          // If basket from queue is null, exit
          if (!basket)
            break;
        }
      }
    }
    // Check if there is any idle propagator in the run manager
    if (!firstTime && multiPropagator) {
      idle = propagator->fRunMgr->GetIdlePropagator();
      if (idle) {
        // Try to steal some work from the run manager
        runmgr->ProvideWorkTo(idle);
      }
    }
    // Start transporting the basket
    MaybeCleanupBaskets(td,basket);
    ++counter;
    ntotransport = basket->GetNinput(); // all tracks to be transported
                                        //      ninput = ntotransport;
    GeantTrack_v &input = basket->GetInputTracks();
    GeantTrack_v &output = *td->fTransported;
    if (!ntotransport)
      goto finish; // input list empty
    //      Geant::Print("","======= BASKET %p with %d tracks counter=%d =======", basket, ntotransport,
    //      counter.load());
    //      basket->Print();
    //      Geant::Print("","==========================================");
    //      propagator->fTracksPerBasket[tid] = ntotransport;
    td->fVolume = 0;
    mat = 0;
    if (!basket->IsMixed()) {
      td->fVolume = basket->GetVolume();
#ifdef USE_VECGEOM_NAVIGATOR
      mat = (Material_t *)td->fVolume->GetMaterialPtr();
#else
      mat = td->fVolume->GetMaterial();
#endif
      if (ntotransport < 257)
        nvect[ntotransport] += ntotransport;
    } else {
      nvect[1] += ntotransport;
    }

    // Select the discrete physics process for all particles in the basket
    if (propagator->fConfig->fUsePhysics) {
      propagator->ProposeStep(ntotransport, input, td);
      // Apply msc for charged tracks; NOTE: we do nothing in this method at the moment
      // propagator->ApplyMsc(ntotransport, input, td);
    }

    ncross = 0;
    generation = 0;

    while (ntotransport) {
      // Interrupt condition here. Work stealing could be also implemented here...
      generation++;
      // Use fNsteps track data to detect geometry anomalies
      for (auto itr=0; itr<ntotransport; ++itr) {
        input.fNstepsV[itr]++;
        if ((input.fStatusV[itr] != kKilled) &&
            (input.fNstepsV[itr] > propagator->fConfig->fNstepsKillThr) &&
            input.fBoundaryV[itr] && (input.fSnextV[itr]<1.e-9)) {
          Error("TransportTracks", "track %d seems to be stuck -> killing it after next step", input.fParticleV[itr]);
          Error("TransportTracks", "Transport will continue, but this is a fatal error");
          input.PrintTrack(itr, "stuck");
          input.fStatusV[itr] = kKilled;
        }
      }

      //         Geant::Print("","====== WorkloadManager:");
      //         input.PrintTracks();
      // Propagate all remaining tracks
      if (basket->IsMixed())
        ncross += input.PropagateTracksScalar(td);
      else
        ncross += input.PropagateTracks(td);
      ntotransport = input.GetNtracks();
    }
    // All tracks are now in the output track vector. Possible statuses:
    // kCrossing - particles crossing boundaries
    // kPhysics - particles reaching the point where the discrete physics process
    //            will happen.
    // kExitingSetup - particles exiting the geometry
    // kKilled - particles that could not advance in geometry after several tries

    // Update track time for all output tracks (this should vectorize)
    nout = output.GetNtracks();
    for (auto itr = 0; itr < nout; ++itr) {
      output.fTimeV[itr]    += output.TimeStep(itr, output.fStepV[itr]);
      // update number-of-interaction-legth-left
      output.fNintLenV[itr] -= output.fStepV[itr]/output.fIntLenV[itr];
    }

    // Post-step actions by continuous processes for all particles. There are no
    // new generated particles at this point.
    if (propagator->fConfig->fUsePhysics) {
      nphys = 0;
      nextra_at_rest = 0;
      // count phyics steps here
      for (auto itr = 0; itr < nout; ++itr)
        if (output.fStatusV[itr] == kPhysics)
          nphys++;

#ifdef USE_REAL_PHYSICS
      propagator->GetPhysicsInterface()->AlongStepAction(mat, output.GetNtracks(), output, nextra_at_rest, td);
#else
      propagator->Process()->Eloss(mat, output.GetNtracks(), output, nextra_at_rest, td);
#endif
      //         if (nextra_at_rest) Geant::Print("","Extra particles: %d", nextra_at_rest);
      // Now we may also have particles killed by energy threshold
      // Do post-step actions on remaining particles
      // to do: group particles per process

      // Discrete processes only
      nphys = output.SortByLimitingDiscreteProcess(); // only those that kPhysics and not continous limit
      if (nphys) {
        // reset number-of-interaction-legth-left
        for (auto itr = 0; itr < nphys; ++itr)
          output.fNintLenV[itr] = -1.0;

        // Do post step actions for particles suffering a given process.
        // Surviving particles are added to the output array
#ifdef USE_REAL_PHYSICS
        propagator->GetPhysicsInterface()->PostStepAction(mat, nphys, output, ntotnext, td);
#else
        // first: sample target and type of interaction for each primary tracks
        propagator->Process()->PostStepTypeOfIntrActSampling(mat, nphys, output, td);

//
// TODO: vectorized final state sampling can be inserted here through
//       a proper interface code that will:
//         1. select the necessary primary tracks based on the alrady
//            sampled interaction type
//         2. take all the member of the selected primary tracks that
//            necessary for sampling the final states
//         3. call the appropriate vector physics code that will
//            perform the physics interaction itself
//         4. update properties of the selected primary tracks,
//            insert secondary tracks into the track vector and set
//            'ntotnext' to the value of the number of secondary tracks
//            inserted to the track vector
//
  #ifdef USE_VECPHYS
        propagator->fVectorPhysicsProcess->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);
  #endif
        // second: sample final states (based on the inf. regarding sampled
        //         target and type of interaction above), insert them into
        //         the track vector, update primary tracks;
        propagator->Process()->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);
#endif
      }
    }
    if (propagator->fStdApplication)
      propagator->fStdApplication->StepManager(output.GetNtracks(), output, td);
    propagator->fApplication->StepManager(output.GetNtracks(), output, td);

    // WP
    #ifdef USE_ROOT
    if (concurrentWrite) {
      while (!(myhitFactory->fOutputsArray[tid].empty())) {
        data = myhitFactory->fOutputsArray[tid].back();
        myhitFactory->fOutputsArray[tid].pop_back();

        tree->Fill();
        // now we can recycle data memory
        myhitFactory->Recycle(data, tid);
      }
      if (tree->GetEntries() > treeSizeWriteThreshold) file->Write();
    }
     #endif

    // Update geometry path for crossing tracks
    ntotnext = output.GetNtracks();

// Normalize directions (should be not needed)
//    for (auto itr=0; itr<ntotnext; ++itr)
//      output.Normalize(itr);

    for (auto itr = 0; itr < ntotnext; ++itr) {
      //output.fNstepsV[itr]++;
      if (output.fStatusV[itr] == kBoundary)
        *output.fPathV[itr] = *output.fNextpathV[itr];
    }
  finish:
    firstTime = false;
    // Check if there are enough transported tracks staying in the same volume
    // to be reused without re-basketizing
    int nreusable = sch->ReusableTracks(output);
    bool reusable = (basket->IsMixed())? false : (nreusable>=propagator->fConfig->fNminReuse);
//    reusable = false;

    if (reusable)
      sch->CopyReusableTracks(output, input, basket->GetThreshold());
    // Remaining tracks need to be re-basketized
    sch->AddTracks(output, ntot, nnew, nkilled, td);
    // Make sure the basket is not recycled before gets released by basketizer
    // This should not happen very often, just when some threads are highly
    // demoted ant the basket makes it through the whole cycle before being fully released
    if (!reusable) {
      // Recycle basket, otherwise keep it for the next iteration
      while (basket->fNused.load())
        ;
      if (basket != bserv) basket->Recycle(td);
      basket = 0; // signal reusability by non-null pointer
    }
    // Update boundary crossing counter
    td->fNcross += ncross;
  }

  // WP
  #ifdef USE_ROOT
  if(concurrentWrite)
    {
      file->Write();
    }
   #endif
  wm->DoneQueue()->push_force(nullptr);
  delete prioritizer;
  // Final reduction of counters
  propagator->fNsteps += td->fNsteps;
  propagator->fNsnext += td->fNsnext;
  propagator->fNphys += td->fNphys;
  propagator->fNmag += td->fNmag;
  propagator->fNsmall += td->fNsmall;
  propagator->fNcross += td->fNcross;

  // If transport is completed, make send the signal to the run manager
  if (runmgr->TransportCompleted())
    runmgr->StopTransport();
  runmgr->GetTDManager()->ReleaseTaskData(td);
  Geant::Print("","=== Thread %d: exiting ===", tid);

  #ifdef USE_ROOT
  if (wm->IsStopped()) wm->MergingServer()->Finish();

  if (concurrentWrite) {
    delete file;
  }

  #endif
  return 0;
}

//______________________________________________________________________________
void *WorkloadManager::TransportTracksCoprocessor(GeantPropagator *prop,TaskBroker *broker) {
  // Thread propagating all tracks from a basket.
  //      char slist[256];
  //      TString sslist;
  //   const int max_idle = 1;
  //   int indmin, indmax;
  static std::atomic<int> counter(0);
  // int ncross;
  int nbaskets;
  int ntotransport;
  // int nextra_at_rest = 0;
  int generation = 0;
  //   int ninjected = 0;
  int nnew = 0;
  int ntot = 0;
  int nkilled = 0;
  int ngcoll = 0;
  unsigned int error = 0;
  GeantBasket *basket = 0;

  GeantPropagator *propagator = prop;
  int basket_size = propagator->fConfig->fNperBasket;
  GeantRunManager *runmgr = propagator->fRunMgr;
  Geant::GeantTaskData *td = runmgr->GetTDManager()->GetTaskData();
  int tid = td->fTid;
  td->AttachPropagator(prop, 0);
  Geant::Print("","=== Worker thread %d created for Coprocessor ===", tid);

  int nworkers = propagator->fNthreads;
  WorkloadManager *wm = propagator->fWMgr;
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  GeantScheduler *sch = wm->GetScheduler();
  int *nvect = sch->GetNvect();
  GeantBasketMgr *prioritizer = new GeantBasketMgr(prop,sch, 0, 0, true);
  td->fBmgr = prioritizer;
  prioritizer->SetThreshold(basket_size);
  prioritizer->SetFeederQueue(feederQ);

  GeantEventServer *evserv = runmgr->GetEventServer();
  int bindex = evserv->GetBindex();
  GeantBasket *bserv = sch->GetBasketManagers()[bindex]->GetNextBasket(td);
  bserv->SetThreshold(basket_size);
  td->fImported = bserv;

#ifndef USE_VECGEOM_NAVIGATOR
  // Create navigator in ROOT case if none serving this thread.
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif

  // Activate events in the server
  //evserv->ActivateEvents();

  // broker->SetPrioritizer(prioritizer);
  while (1) {

    // Call the feeder if in priority mode
    auto feedres = wm->CheckFeederAndExit();
    if (feedres == FeederResult::kStop) {
       break;
    }

    // ::Info("GPU","Waiting (1) for next available stream.");
    // TaskBroker::Stream stream = broker->GetNextStream();
    // if (!stream)
    //   break;

    if (wm->FeederQueue()->empty()) {
      // There is no work to be done for now, let's just run what we have
      if (0 != broker->launchTask()) {
        // We ran something, let wait for the next free stream,
        MaybeCleanupBaskets(td);
        continue;
      } else {
        // We had nothing to run at all .... we need to wait for
        // more data ....
      }
    }
    nbaskets = feederQ->size_async();
    if (nbaskets > nworkers)
      ngcoll = 0;
    // If prioritizers have work, just do it
    if ((basket = broker->GetBasketForTransport(*td))) {
      ngcoll = 0;
    } else {
      if (nbaskets < 1 && !evserv->HasTracks() ) {
        sch->GarbageCollect(td);
        ngcoll++;
      }
      // Too many garbage collections - enter priority mode
      if ((ngcoll > 5) && (propagator->GetNworking() <= 1)) {
        ngcoll = 0;
        while (!sch->GarbageCollect(td, true) &&
               !feederQ->size_async() &&
               !prioritizer->HasTracks() &&
               !evserv->HasTracks())
          ;
      }
      if (wm->FeederQueue()->try_pop(basket)) {
         // We got a basket, let's go.
      } else {
         // No basket in the queue, let's flush our oustanding work
         if (0 != broker->launchTask(/* wait= */ true)) {
            // We ran something, new basket might be available.
            continue;
         } else {
           // Get basket from the generator
           if (evserv->HasTracks()) {
           // In the initial phase we distribute a fair share of baskets to all propagators
             if (!evserv->IsInitialPhase() ||
                propagator->fNbfeed < runmgr->GetInitialShare()) {
               ntotransport = evserv->FillBasket(bserv->GetInputTracks(), basket_size, error);
               if (ntotransport) basket = bserv;
             }
           }
           // We have nothing, so let's wait.

           if (propagator->fCompleted || wm->fStopped)
             break;
           propagator->fNidle++;
           wm->FeederQueue()->wait_and_pop(basket);
           propagator->fNidle--;
         }
      }
    }
    // Check exit condition: null basket in the queue
    if (!basket) {
      // All work should have been done, the broker's queue should
      // be empty.
      // assert( worker queue empty );
      break;
    }
    MaybeCleanupBaskets(td,basket);
    ++counter;

    // if (!stream) {
    //   ::Info("GPU", "Waiting (2) for next available stream.");
    //   stream = broker->GetNextStream();
    // }
    // lastToClear = false;
    // if (!basket) {
    //   if (0 != broker->launchTask(/* wait= */ true)) {
    //     // We ran something, new basket might be available.
    //     continue;
    //   } else {
    //     break;
    //   }
    // }

    ntotransport = basket->GetNinput(); // all tracks to be transported
    // ninput = ntotransport;
    //    GeantTrack_v &input = basket->GetInputTracks();
    GeantTrack_v &output = *td->fTransported;
    if (!ntotransport)
      goto finish; // input list empty
    //      Geant::Print("","======= BASKET %p with %d tracks counter=%d =======", basket, ntotransport,
    //      counter);
    //      basket->Print();
    //      Geant::Print("","==========================================");
    //      propagator->fTracksPerBasket[tid] = ntotransport;
    td->fVolume = 0;
    // mat = 0;
    if (!basket->IsMixed()) {
      td->fVolume = basket->GetVolume();
#ifdef USE_VECGEOM_NAVIGATOR
      // mat = (Material_t *)td->fVolume->GetMaterialPtr();
#else
      // mat = td->fVolume->GetMaterial();
#endif
      if (ntotransport < 257)
        nvect[ntotransport] += ntotransport;
    } else {
      nvect[1] += ntotransport;
    }

    // Record tracks
    // ninput = ntotransport;
    //      if (counter==1) input.PrintTracks();
    for (int itr = 0; itr < ntotransport; itr++) {
      // iev[itr] = input.fEventV[itr];
      // itrack[itr] = input.fParticleV[itr];
      // crt[itr] = input.fPathV[itr];
      // nxt[itr] = input.fNextpathV[itr];
      //if (isnan(input.fXdirV[itr])) {
      //  Geant::Error("Before transport", "Error: track %d has NaN", itr);
      //}
    }

    // Select the discrete physics process for all particles in the basket
    // Moved to CoprocessorBroker kernel.
    // if (propagator->fConfig->fUsePhysics) {
    //   propagator->ProposeStep(ntotransport, input, td);
    //   // Apply msc for charged tracks
    //   propagator->ApplyMsc(ntotransport, input, td);
    // }

    // ncross = 0;
    generation = 0;

    while (ntotransport) {
      // Interrupt condition here. Work stealing could be also implemented here...
      generation++;
      // Propagate all remaining tracks
      // NOTE: need to deal with propagator->fConfig->fUsePhysics
      broker->runTask(*td, *basket); // ntotransport, basket_sch->GetNumber(), propagator->fTracks, particles);
      ntotransport = 0;
      // ncross += input.PropagateTracks(output);
      // ntotransport = input.GetNtracks();
    }
    // All tracks are now in the output track vector. Possible statuses:
    // kCrossing - particles crossing boundaries
    // kPhysics - particles reaching the point where the discrete physics process
    //            will happen.
    // kExitingSetup - particles exiting the geometry
    // kKilled - particles that could not advance in geometry after several tries

    {
       auto noutput = basket->GetNinput();
       for (int itr = 0; itr < noutput; itr++) {
       #ifdef USE_ROOT
          if (TMath::IsNaN(output.fXdirV[itr])) {
       #else
          if (std::isnan(output.fXdirV[itr])) {
       #endif
             Geant::Error("TransportTracksCoprocessor","Track %d has NaN", itr);
          }
       }
    }
    // Note: Need to apply the code if (propagator->fConfig->fUsePhysics)
    // to the tracks after they came back.

  finish:
//      basket->Clear();
//      Geant::Print("","======= BASKET(tid=%d): in=%d out=%d =======", tid, ninput, basket->GetNoutput());
    /* int ninjected = */ sch->AddTracks(output, ntot, nnew, nkilled, td);
    if (prioritizer->HasTracks()) {
       // We know that all the left over tracks injected here can not be handle
       // by the coprocessor, so we much send it for somebody else to process:

       GeantBasket *cbasket = prioritizer->GetCBasket();
       assert(cbasket->TryDispatch());
       prioritizer->SetCBasket(prioritizer->GetNextBasket(td));
       prioritizer->Push(cbasket, /*priority=*/ false, td);
       // prioritizer->BookBasket(td);
       // prioritizer->GarbageCollect(td);
    }
    //      Geant::Print("","thread %d: injected %d baskets", tid, ninjected);
    // wm->TransportedQueue()->push(basket);
    (void)ntot;
    (void)nnew;
    (void)nkilled;
//    sched_locker.StartOne();
    // Make sure the basket is not recycled before gets released by basketizer
    // This should not happen vey often, just when some threads are highly
    // demoted ant the basket makes it through the whole cycle before being fully released
    while (basket->fNused.load()) ;
    if (basket != bserv) basket->Recycle(td);
  }
  // delete prioritizer;
  wm->DoneQueue()->push_force(nullptr);
  delete prioritizer;
  // Final reduction of counters
  propagator->fNsteps += td->fNsteps;
  propagator->fNsnext += td->fNsnext;
  propagator->fNphys += td->fNphys;
  propagator->fNmag += td->fNmag;
  propagator->fNsmall += td->fNsmall;
  runmgr->GetTDManager()->ReleaseTaskData(td);
  Geant::Print("","=== Coprocessor Thread %d: exiting === Processed %ld", tid, broker->GetTotalWork());
  return 0;
}

//______________________________________________________________________________
void *WorkloadManager::GarbageCollectorThread(GeantPropagator *prop) {
#ifdef USE_ROOT
  // This threads can be triggered to do garbage collection of unused baskets
  static double rsslast = 0;
  double rss;
  ProcInfo_t procInfo;
  const double MByte = 1024.;
  const double_t thr_increase = 1.05;
  GeantPropagator *propagator = prop;
  WorkloadManager *wm = propagator->fWMgr;
  int nthreads = propagator->fNthreads;
  double threshold = propagator->fConfig->fMaxRes;
  double virtlimit = propagator->fConfig->fMaxVirt;
  if (threshold == 0 && virtlimit == 0)
    return 0;
  //  condition_locker &gbc_locker = wm->GetGbcLocker();
  while (1) {
    //    gbc_locker.Wait();
    if (wm->IsStopped())
      break;
    // todo: cleaning with thread based recycle queues
    gSystem->GetProcInfo(&procInfo);
    rss = procInfo.fMemResident / MByte;
    double vsize = procInfo.fMemVirtual / MByte;
    Geant::Print("GC thread","### Checking mem: %g MBytes", rss);
    bool needClean = false;
    if (threshold && rss > threshold && (rss / rsslast > thr_increase)) {
      rsslast = rss;
      needClean = true;
    } else if (virtlimit && vsize > virtlimit) {
      needClean = true;
    }
    if (needClean) {
      for (int tid = 0; tid < nthreads; tid++) {
        GeantTaskData *td = propagator->fRunMgr->GetTDManager()->GetTaskData(tid);
        td->SetToClean(true);
      }
    }
    gSystem->Sleep(1000); // millisec
  }
#else
  (void)prop;
#endif
  return 0;
}



//______________________________________________________________________________
void *WorkloadManager::MonitoringThread(GeantPropagator* prop) {
#ifdef USE_ROOT
  // Thread providing basic monitoring for the scheduler.
  const double MByte = 1024.;
  Geant::Info("MonitoringThread","Started monitoring ...");
  std::cout << "isBatch = " << gROOT->IsBatch() << std::endl;
  GeantPropagator *propagator = prop;
  WorkloadManager *wm = propagator->fWMgr;
  int nmon = prop->fConfig->GetMonFeatures();
  if (!nmon)
    return 0;
  int dmon = 0.5 * nmon + nmon % 2;
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  int ntotransport;
  ProcInfo_t procInfo;
  double rss, rssmax = 0;
  double nmem[101] = {0};
  GeantScheduler *sch = wm->GetScheduler();
  int nvol = sch->GetNvolumes();
  int *nvect = sch->GetNvect();
  int nthreads = wm->GetNthreads();
  int *nworking = new int[nthreads + 1];
  memset(nworking, 0, (nthreads + 1) * sizeof(int));
  GeantBasketMgr **bmgr = sch->GetBasketManagers();
  TCanvas *cmon = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("cscheduler");
  if (nmon == 1)
    cmon->Divide(1, 1);
  else
    cmon->Divide(2, dmon);
  TH1I *hqueue = 0;
  int nqueue[101] = {0};
  int ipad = 0;
  if (propagator->fConfig->IsMonitored(GeantConfig::kMonQueue)) {
    hqueue = new TH1I("hqueue", "Work queue load", 100, 0, 100);
    hqueue->SetFillColor(kRed);
    hqueue->SetLineColor(0);
    hqueue->SetStats(false);
    cmon->cd(++ipad);
    hqueue->Draw();
  }
  TH1F *hmem = 0;
  if (propagator->fConfig->IsMonitored(GeantConfig::kMonMemory)) {
    hmem = new TH1F("hmem", "Resident memory [MB]", 100, 0, 100);
    // hmem->SetFillColor(kMagenta);
    hmem->SetLineColor(kMagenta);
    hmem->SetStats(false);
    cmon->cd(++ipad);
    hmem->Draw();
  }
  TH1I *hbaskets = 0;
  TH1I *hbused = 0;
  if (propagator->fConfig->IsMonitored(GeantConfig::kMonBasketsPerVol)) {
    hbaskets = new TH1I("hbaskets", "Baskets per volume", nvol, 0, nvol);
    hbaskets->SetFillColor(kBlue);
    hbaskets->SetLineColor(0);
    hbaskets->SetStats(false);
    //    hbaskets->GetYaxis()->SetRangeUser(1,10000);
    hbused = new TH1I("hbused", "Baskets per volume", nvol, 0, nvol);
    hbused->SetFillColor(kRed);
    hbused->SetFillStyle(3001);
    hbused->SetLineColor(0);
    hbused->SetStats(false);
    cmon->cd(++ipad)->SetLogy();
    hbaskets->Draw();
    hbused->Draw("SAME");
  }
  TH1I *hvectors = 0;
  if (propagator->fConfig->IsMonitored(GeantConfig::kMonVectors)) {
    hvectors = new TH1I("hvectors", "Tracks in vectors of given size", 257, 0, 257);
    hvectors->SetFillColor(kBlue);
    hvectors->SetLineColor(0);
    hvectors->SetStats(false);
    cmon->cd(++ipad)->SetLogy();
    hvectors->Draw();
  }
  TH1I *htracksmax = 0;
  TH1I *htracks = 0;
  int nbuffered = propagator->fNbuff;
  if (propagator->fConfig->IsMonitored(GeantConfig::kMonTracksPerEvent)) {
    htracksmax = new TH1I("htracksmax", "Tracks in flight", nbuffered, 0, nbuffered);
    htracksmax->SetFillColor(kBlue);
    htracksmax->SetLineColor(0);
    htracksmax->SetStats(false);
    htracks = new TH1I("htracks", "Tracks in flight", nbuffered, 0, nbuffered);
    htracks->SetFillColor(kRed);
    htracks->SetFillStyle(3001);
    htracks->SetLineColor(0);
    htracks->SetStats(false);
    cmon->cd(++ipad)->SetLogy();
    htracksmax->Draw();
    htracks->Draw("SAME");
  }
  TH1I *htrackstot = 0;
  int ntrackstot[101] = {0};
  if (propagator->fConfig->IsMonitored(GeantConfig::kMonTracks)) {
    htrackstot = new TH1I("htrackstot", "Total number of tracks alive", 100, 0, 100);
    htrackstot->SetFillColor(kRed);
    htrackstot->SetLineColor(0);
    htrackstot->SetStats(false);
    cmon->cd(++ipad);
    htrackstot->Draw();
  }
  cmon->Update();
  double stamp = 0.;
  int i, j, bin;
  int nmaxtot;
  while (1) { // exit condition here
    i = int(stamp);
    ipad = 0;
    gSystem->Sleep(50); // millisec
    // Fill histograms
    if (stamp > 100) {
      if (hqueue) {
        ntotransport = feederQ->size_async();
        memmove(nqueue, &nqueue[1], 99 * sizeof(int));
        nqueue[99] = ntotransport;
        hqueue->GetXaxis()->Set(100, stamp - 100, stamp);
        for (j = 0; j < 100; j++)
          hqueue->SetBinContent(j + 1, nqueue[j]);
      }
      if (hmem) {
        gSystem->GetProcInfo(&procInfo);
        rss = procInfo.fMemResident / MByte;
        if (rss > rssmax) {
          rssmax = rss;
          printf("RSS = %g MB\n", rssmax);
        }
        memmove(nmem, &nmem[1], 99 * sizeof(double));
        nmem[99] = rss;
        hmem->GetXaxis()->Set(100, stamp - 100, stamp);
        for (j = 0; j < 100; j++)
          hmem->SetBinContent(j + 1, nmem[j]);
      }
      if (htrackstot) {
        // Count tracks for all event slots
        int ntr = 0;
        for (int slot = 0; slot < propagator->fNbuff; slot++)
          ntr += propagator->fEvents[slot]->GetNinflight();
        memmove(ntrackstot, &ntrackstot[1], 99 * sizeof(int));
        ntrackstot[99] = ntr;
        htrackstot->GetXaxis()->Set(100, stamp - 100, stamp);
        for (j = 0; j < 100; j++)
          htrackstot->SetBinContent(j + 1, ntrackstot[j]);
      }
    } else {
      if (hqueue) {
        ntotransport = feederQ->size_async();
        nqueue[i] = ntotransport;
        hqueue->SetBinContent(i + 1, nqueue[i]);
      }
      if (hmem) {
        gSystem->GetProcInfo(&procInfo);
        rss = procInfo.fMemResident / MByte;
        if (rss > rssmax) {
          rssmax = rss;
          printf("RSS = %g MB\n", rssmax);
        }
        nmem[i] = rss;
        hmem->SetBinContent(i + 1, nmem[i]);
      }
      if (htrackstot) {
        // Count tracks for all event slots
        int ntr = 0;
        for (int slot = 0; slot < propagator->fNbuff; slot++)
          ntr += propagator->fEvents[slot]->GetNinflight();
        ntrackstot[i] = ntr;
        htrackstot->SetBinContent(i + 1, ntrackstot[i]);
      }
    }
    if (hbaskets) {
      for (j = 0; j < nvol; j++) {
        bin = hbaskets->GetXaxis()->FindBin(bmgr[j]->GetVolume()->GetName());
        hbaskets->SetBinContent(bin, bmgr[j]->GetNbaskets());
        hbused->SetBinContent(bin, bmgr[j]->GetNused());
      }
      int nbaskets_mixed = 0;
      int nused_mixed = 0;
      for (j = 0; j < nthreads; j++) {
        GeantTaskData *td = propagator->fRunMgr->GetTDManager()->GetTaskData(j);
        nbaskets_mixed += td->fBmgr->GetNbaskets();
        nused_mixed += td->fBmgr->GetNused();
      }
      bin = hbaskets->GetXaxis()->FindBin("MIXED");
      hbaskets->SetBinContent(bin, nbaskets_mixed);
      hbused->SetBinContent(bin, nused_mixed);
    }
    if (hvectors) {
      int vmax = 0;
      int ymax = 0;
      // mixed tracks
      // hvectors->SetBinContent(1, nvect[1]);
      for (j = 1; j < 257; ++j) {
        if (nvect[j] > 0) {
          if (j > vmax)
            vmax = j;
        }
        if (nvect[j] > ymax)
          ymax = nvect[j];
        hvectors->SetBinContent(j + 1, nvect[j]);
      }
      hvectors->GetXaxis()->SetRangeUser(-5, vmax);
      hvectors->GetYaxis()->SetRangeUser(1, 1.2 * ymax);
    }
    if (htracksmax) {
      nmaxtot = 1;
      for (j = 0; j < nbuffered; j++) {
        GeantEvent *evt = propagator->fEvents[j];
        int nmax = evt->GetNmax();
        nmaxtot = max<int>(nmax, nmaxtot);
        htracksmax->SetBinContent(j + 1, nmax);
        htracks->SetBinContent(j + 1, evt->GetNinflight());
      }
      htracksmax->GetYaxis()->SetRangeUser(1, 10 * nmaxtot);
    }
    // Now draw all histograms
    if (hqueue) {
      cmon->cd(++ipad);
      hqueue->Draw();
    }
    if (hmem) {
      cmon->cd(++ipad);
      hmem->Draw();
    }
    if (hbaskets) {
      cmon->cd(++ipad);
      hbaskets->LabelsOption(">");
      hbaskets->GetXaxis()->SetRangeUser(0, 10);
      hbaskets->Draw();
      hbused->Draw("SAME");
    }
    if (hvectors) {
      cmon->cd(++ipad);
      hvectors->Draw();
    }
    if (htracksmax) {
      cmon->cd(++ipad);
      htracksmax->GetXaxis()->SetTitle("Event slot");
      htracksmax->Draw();
      htracks->Draw("SAME");
    }
    if (htrackstot) {
      cmon->cd(++ipad);
      htrackstot->Draw();
    }
    cmon->Modified();
    cmon->Update();
    stamp += 1;
    if (wm->IsStopped())
      break;
  }
  delete[] nworking;
  double sumvect = 0.;
  for (i = 0; i < 257; ++i)
    sumvect += nvect[i];
  Geant::Info("Monitoring stopped","=== Percent of tracks transported in single track mode: %g%%", 100. * nvect[1] / sumvect);
  // Sleep a bit to let the graphics finish
  gSystem->Sleep(100); // millisec
#else
  (void)prop;
#endif
  return 0;
}

//______________________________________________________________________________
void *WorkloadManager::OutputThread(GeantPropagator* prop) {
  // Thread providing basic output for the scheduler.

  Geant::Info("OutputThread","%s","=== Output thread created ===");
  #ifdef USE_ROOT

  if (prop->fConfig->fConcurrentWrite) {
    Printf(">>> Writing concurrently to MemoryFiles");

    prop->fWMgr->MergingServer()->Listen();
  }
  else {
    TFile file("hits.root", "RECREATE");
    WorkloadManager *wm = prop->fWMgr;

    GeantBlock <MyHit> *data = 0;

    TTree *tree = new TTree("Hits", "Simulation output");

    tree->Branch("hitblocks", "GeantBlock<MyHit>", &data);

    GeantFactoryStore *factoryStore = GeantFactoryStore::Instance();
    GeantFactory <MyHit> *myhitFactory = factoryStore->GetFactory<MyHit>(16, prop->fRunMgr->GetNthreadsTotal());


    while (!(wm->IsStopped()) || myhitFactory->fOutputs.size() > 0) {
      // Geant::Print("","size of queue from output thread %zu", myhitFactory->fOutputs.size());

      if (myhitFactory->fOutputs.size() > 0) {
        while (myhitFactory->fOutputs.try_pop(data)) {
          // myhitFactory->fOutputs.wait_and_pop(data);
          // Geant::Print("","Popping from queue of %zu", myhitFactory->fOutputs.size() + 1);
          // if(data) std::cout << "size of the block in the queue " << data->Size() << std::endl;

          tree->Fill();
          myhitFactory->Recycle(data);
        }
      }
    }
    tree->Write();
    file.Close();
  }

  Geant::Info("OutputThread", "=== Output thread finished ===");

  return 0;

  #else
    (void)prop;
    printf("%s","=== ROOT is disabled - output thread did nothing ===");
  #endif
    return 0;
}

#ifdef USE_ROOT
//______________________________________________________________________________
void *WorkloadManager::StartROOTApplication()
{
  TApplication *app = new TApplication("GeantV monitoring", NULL, NULL);
  app->Run();
  return 0;
}
#endif

} // GEANT_IMPL_NAMESPACE
} // Geant
