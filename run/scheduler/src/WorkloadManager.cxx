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
#include "GeantTaskData.h"
#include "PhysicsInterface.h"
#include "PhysicsProcessOld.h"
#include "GeantEvent.h"
#include "EventSet.h"
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
WorkloadManager *WorkloadManager::NewInstance(GeantPropagator *prop, int nthreads) {
  // Return singleton instance.
  if (!nthreads || !prop) {
    Geant::Error("WorkloadManager::NewInstance", "%s", "You should provide number of threads.");
    return 0;
  }
  return new WorkloadManager(nthreads,prop);
}

void WorkloadManager::SetTaskBroker(TaskBroker *broker) {
  // Register the broker if it is valid.
  if (broker && broker->IsValid())
    fBroker = broker;
}

//______________________________________________________________________________
bool WorkloadManager::StartTasks() {
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
    // fListThreads.emplace_back(WorkloadManager::TransportTracksCoprocessor, prop, fBroker);
    ith += fBroker->GetNstream() + 1;
    if (ith == fNthreads && fBroker->IsSelective()) {
       Geant::Fatal("WorkloadManager::StartThreads","All %d threads are used by the coprocessor broker but it can only process a subset of particles.",fNthreads);
       return false;
    }
  }

  // Start CPU transport threads (static mode)
  for (; ith < fNthreads; ith++)
    fListThreads.emplace_back(WorkloadManager::TransportTracksV3, prop);

  return true;
}

//______________________________________________________________________________
void WorkloadManager::WaitWorkers() {
  // Waiting point for the main thread until work gets done.
  int ntowait = fNthreads;
  if (fBroker)
    ntowait -= fBroker->GetNstream();
  Basket *signal;
  while (ntowait) {
    fDoneQ->wait_and_pop(signal);
    ntowait--;
    Geant::Print("", "=== %d workers finished", fNthreads - ntowait);
  }
  //   fBasketGeneration++;
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
/** @brief Call Feeder (if needed) and check exit condition. */
WorkloadManager::FeederResult WorkloadManager::CheckFeederAndExit() {

//  if (!prioritizer.HasTracks() && (propagator.fRunMgr->GetNpriority() || propagator.GetNworking() == 1)) {
//    bool didFeeder = propagator.fRunMgr->Feeder(&td) > 0;

  // Check exit condition
  if (fPropagator->fRunMgr->TransportCompleted()) {
    Stop();
    return FeederResult::kStop;
  }
  return FeederResult::kNone;
}

//______________________________________________________________________________
bool WorkloadManager::TransportTracksTask(EventSet *workload, GeantTaskData *td) {
// Re-entrant main transport method. This will co-operate with other identical
// concurrent tasks (basketizing + transporting tracks for all events available
// in the event server). The method will exit if it triggered finishing any of
// the event sets registered in the event server.
// Remarks:
//   - this is a task mode, NUMA not enabled
//   - the task data has to be pre-booked using RunManager::BookTransportTask
   if (!td) {
    Error("TransportTracksTask", "No task data object available!!!");
    return false;
  }

  GeantPropagator *propagator = td->fPropagator;
  GeantRunManager *runmgr = propagator->fRunMgr;
#ifndef USE_VECGEOM_NAVIGATOR
  // If we use ROOT make sure we have a navigator here
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif
/*** STEPPING LOOP ***/
  bool flush = false;
  while (!workload->IsDone()) {
    // Check if simulation has finished
    auto feedres = PreloadTracksForStep(td);
    flush = (feedres == FeederResult::kNone) ||
            (feedres == FeederResult::kError) ||
            workload->IsDone();
    SteppingLoop(td, flush);
    if (workload->IsDone()) break;
    if (flush) {
      // There is no more work in the server but the workload is not yet
      // completed. We cannot return as the workload is not done, so we need to
      // put the thread on hold. This should not be significant overhead given
      // that it can only happen if the worker has finished its own work and there
      // is no other work from the server.
      
      workload->SleepUntilDone();
      assert(workload->IsDone());
      continue;
    }
  }
  
  delete workload;

  // Release the task data
  runmgr->GetTDManager()->ReleaseTaskData(td);
  return true;
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
      ninjected = evserv->FillStackBuffer(td->fStackBuffer, td->fPropagator->fConfig->fNperBasket, td, error);
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
// The main stepping loop over simulation stages. Called in flush mode when no tracks are
// available from the event server. The flush mode works as following: if input tracks are
// not available for a stage, the staged is processed in flush mode, otherwise the input
// tracks are processed normally. The flush mode stops if tracks are available again in the
// server.

  int nprocessed = 0;
  int nproc = 0;
  int ninput = 0;
  int istage = 0;
  int ninjected = FlushOneLane(td);
  // Flush the input buffer lane by lane, starting with higher generations.
  // Exit loop when buffer is empty or when the stages are cosidered flushed
  while ( ninjected || flush ) {
    while (1) {
      // How many particles at input of current stage?
      SimulationStage *stage = td->fPropagator->fStages[istage];
      int nstart = td->fStageBuffers[istage]->size();
      ninput += nstart;
      // If there is input just process normally
      if (nstart) {
        nproc += stage->Process(td);
      } else {
        // If flush mode requested and the stage is basketized, process in flush mode
        if (flush && stage->IsBasketized())
          nproc += stage->FlushAndProcess(td);
      }
      nprocessed += nproc;
      // Go to next stage
      istage = (istage + 1) % kNstages;
      if (istage == 0) {
        // Checkpoint
        // If no activity in the loop inject next lane from the buffer
        if (ninput == 0 && nproc == 0) {
          ninjected = FlushOneLane(td);
          if (!ninjected) return nprocessed;
          break;
        }
        // TO DO: In flush mode, check regularly the event server
        ninput = 0;
        nproc = 0;
      }
    }
  }
  return nprocessed;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
