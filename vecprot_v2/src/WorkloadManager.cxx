#include "WorkloadManager.h"
#include "Geant/Error.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TBits.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TThread.h"
#include "GeantTrack.h"
#include "GeantBasket.h"
#include "GeantOutput.h"
#include "GeantTaskData.h"
#include "PhysicsProcess.h"
#include "GeantScheduler.h"
#include "GeantEvent.h"
#include "GeantVApplication.h"
#if USE_VECGEOM_NAVIGATOR == 1
#include "management/GeoManager.h"
#include "volumes/Medium.h"
#else
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#endif
#include "TaskBroker.h"

using namespace Geant;

ClassImp(WorkloadManager)

    WorkloadManager *WorkloadManager::fgInstance = 0;

//______________________________________________________________________________
WorkloadManager::WorkloadManager(Int_t nthreads)
    : fNthreads(nthreads), fNbaskets(0), fBasketGeneration(0), fNbasketgen(0), fNidle(nthreads), fNminThreshold(10),
      fNqueued(0), fBtogo(0), fSchId(nthreads), fStarted(false), fStopped(false), fFeederQ(0), fTransportedQ(0),
      fDoneQ(0), fListThreads(0), fFlushed(false), fFilling(false), fMonQueue(0), fMonMemory(0), fMonBasketsPerVol(0),
      fMonVectors(0), fMonConcurrency(0), fMonTracksPerEvent(0), fMonTracks(0), fScheduler(0), fBroker(0), fWaiting(0),
      fSchLocker(), fGbcLocker(), fLastEvent(0) {
  // Private constructor.
  fFeederQ = new Geant::priority_queue<GeantBasket *>(1 << 16);
  fTransportedQ = new Geant::priority_queue<GeantBasket *>(1 << 16);
  fDoneQ = new Geant::priority_queue<GeantBasket *>(1 << 10);
  fgInstance = this;
  fScheduler = new GeantScheduler();
  fWaiting = new Int_t[fNthreads + 1];
  memset(fWaiting, 0, (fNthreads + 1) * sizeof(Int_t));
}

//______________________________________________________________________________
WorkloadManager::~WorkloadManager() {
  // Destructor.
  delete fFeederQ;
  delete fTransportedQ;
  delete fDoneQ;
  delete fScheduler;
  delete[] fWaiting;
  //   delete fNavStates;
  fgInstance = 0;
}

//______________________________________________________________________________
void WorkloadManager::CreateBaskets() {
  // Create the array of baskets
  VolumePath_t *blueprint = 0;
  Int_t maxdepth = GeantPropagator::Instance()->fMaxDepth;
  Printf("Max depth: %d", maxdepth);
  blueprint = VolumePath_t::MakeInstance(maxdepth);
  //   fNavStates = new GeantObjectPool<VolumePath_t>(1000*fNthreads, blueprint);
  //   fNavStates = new rr_pool<VolumePath_t>(16*fNthreads, 1000, blueprint);
  fScheduler->CreateBaskets();
  VolumePath_t::ReleaseInstance(blueprint);

  if (fBroker)
    fBroker->CreateBaskets();
}

//______________________________________________________________________________
WorkloadManager *WorkloadManager::Instance(Int_t nthreads) {
  // Return singleton instance.
  if (fgInstance)
    return fgInstance;
  if (!nthreads) {
    ::Error("WorkloadManager::Instance", "No instance yet so you should provide number of threads.");
    return 0;
  }
  return new WorkloadManager(nthreads);
}

//______________________________________________________________________________
void WorkloadManager::Print(Option_t *) const {
  //
}

void WorkloadManager::SetTaskBroker(TaskBroker *broker) {
  // Register the broker if it is valid.
  if (broker && broker->IsValid())
    fBroker = broker;
}

#if USE_VECGEOM_NAVIGATOR == 1
//______________________________________________________________________________
Bool_t WorkloadManager::LoadGeometry(vecgeom::VPlacedVolume const *const volume) {
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
void WorkloadManager::StartThreads() {
  // Start the threads
  fStarted = true;
  if (fListThreads)
    return;
  fListThreads = new TList();
  fListThreads->SetOwner();
  Int_t ith = 0;
  TThread *t;
  if (fBroker) {
    if (fBroker->GetNstream() > fNthreads) {
      ::Fatal("StartThreads", "The task broker is using too many threads (%d out of %d)", fBroker->GetNstream(),
              fNthreads);
    }
    Printf("Running with a coprocessor broker.");
    t = new TThread(WorkloadManager::TransportTracksCoprocessor, fBroker);
    fListThreads->Add(t);
    t->Run();
    ith += fBroker->GetNstream() + 1;
  }
  // Start CPU transport threads
  for (; ith < fNthreads; ith++) {
    t = new TThread(WorkloadManager::TransportTracks);
    fListThreads->Add(t);
    t->Run();
  }
  //   gSystem->Sleep(1000);
  // Start scheduler(s)
  //  t = new TThread(WorkloadManager::MainScheduler);
  //  fListThreads->Add(t);
  //  t->Run();
  // Start monitoring thread
  if (GeantPropagator::Instance()->fUseMonitoring) {
    t = new TThread(WorkloadManager::MonitoringThread);
    fListThreads->Add(t);
    t->Run();
  }
  // Start garbage collector
  if (GeantPropagator::Instance()->fMaxRes > 0) {
    t = new TThread(WorkloadManager::GarbageCollectorThread);
    fListThreads->Add(t);
    t->Run();
  }
}

//______________________________________________________________________________
void WorkloadManager::JoinThreads() {
  //
  int tojoin = fNthreads;
  if (fBroker)
    tojoin -= fBroker->GetNstream();
  for (Int_t ith = 0; ith < tojoin; ith++)
    fFeederQ->push(0);
  for (Int_t ith = 0; ith < tojoin; ith++)
    ((TThread *)fListThreads->At(ith))->Join();
  // Join scheduler
  //  ((TThread *)fListThreads->At(tojoin))->Join();
  // Join monitoring thread
  if (GeantPropagator::Instance()->fUseMonitoring) {
    ((TThread *)fListThreads->At(tojoin++))->Join();
  }
  // Join garbage collector
  if (GeantPropagator::Instance()->fMaxRes > 0)
    ((TThread *)fListThreads->At(tojoin++))->Join();
}

//______________________________________________________________________________
void WorkloadManager::WaitWorkers() {
  // Waiting point for the main thread until work gets done.
  fFilling = false;
  Int_t ntowait = fNthreads;
  if (fBroker)
    ntowait -= fBroker->GetNstream();
  GeantBasket *signal;
  while (ntowait) {
    fDoneQ->wait_and_pop(signal);
    ntowait--;
    Printf("=== %d workers finished", fNthreads - ntowait);
  }
  //   fBasketGeneration++;
}

//______________________________________________________________________________
void *WorkloadManager::MainScheduler(void *) {
  // Garbage collector thread, called by a single thread.
  Printf("=== Scheduler: stopping threads and exiting ===");
  return 0;
}

//______________________________________________________________________________
void *WorkloadManager::TransportTracks(void *) {
  // Thread propagating all tracks from a basket.
  //      char slist[256];
  //      TString sslist;
  //   const Int_t max_idle = 1;
  //   Int_t indmin, indmax;
  static std::atomic<int> counter(0);
  Int_t ntotnext, ncross, nbaskets;
  Int_t ntotransport;
  Int_t nextra_at_rest = 0;
  Int_t generation = 0;
  //   Int_t ninjected = 0;
  Int_t nnew = 0;
  Int_t ntot = 0;
  Int_t nkilled = 0;
  Int_t nphys = 0;
  Int_t ngcoll = 0;
  GeantBasket *basket = 0;
  Int_t tid = Instance()->ThreadId();
  Printf("=== Worker thread %d created ===", tid);
  GeantPropagator *propagator = GeantPropagator::Instance();
  Geant::GeantTaskData *td = propagator->fThreadData[tid];
  td->fTid = tid;
  Int_t nworkers = propagator->fNthreads;
  WorkloadManager *wm = WorkloadManager::Instance();
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  GeantScheduler *sch = wm->GetScheduler();
  Int_t *nvect = sch->GetNvect();
  GeantBasketMgr *prioritizer = new GeantBasketMgr(sch, 0, 0, true);
  td->fBmgr = prioritizer;
  prioritizer->SetThreshold(propagator->fNperBasket);
  prioritizer->SetFeederQueue(feederQ);
  // Start the feeder
  propagator->Feeder(td);
  TGeoMaterial *mat = 0;
  Int_t *waiting = wm->GetWaiting();
//  condition_locker &sched_locker = wm->GetSchLocker();
//  condition_locker &gbc_locker = wm->GetGbcLocker();
// The last worker wakes the scheduler
//  if (tid == nworkers-1) sched_locker.StartOne();
//   Int_t nprocesses = propagator->fNprocesses;
//   Int_t ninput, noutput;
//   Bool_t useDebug = propagator->fUseDebug;
//   Printf("(%d) WORKER started", tid);
// Create navigator if none serving this thread.
#ifdef USE_VECGEOM_NAVIGATOR
// Suppose I do not need a new navigator for VecGeom... otherwise how it ever worked... ?
#else
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif
  //   Int_t iev[500], itrack[500];
  // TGeoBranchArray *crt[500], *nxt[500];
  while (1) {
    // Call the feeder if in priority mode
    if (!prioritizer->HasTracks() && (propagator->GetNpriority() || wm->GetNworking() == 1)) {
      if (propagator->Feeder(td))
        ngcoll = 0;
      // Check exit condition
      if (propagator->TransportCompleted()) {
        for (Int_t i = 0; i < nworkers; i++)
          wm->FeederQueue()->push(0);
        wm->TransportedQueue()->push(0);
        wm->Stop();
        //         sched_locker.StartOne(); // signal the scheduler who has to exit
        //         gbc_locker.StartOne();
        break;
      }
    }
    waiting[tid] = 1;
    nbaskets = feederQ->size_async();
    if (nbaskets > nworkers)
      ngcoll = 0;
    // If prioritizer has work, just do it
    if (prioritizer->HasTracks()) {
      basket = prioritizer->GetBasketForTransport(td);
      ngcoll = 0;
    } else {
      if ((nbaskets < 1) && (!propagator->IsFeeding())) {
        sch->GarbageCollect(td);
        ngcoll++;
      }
      // Too many garbage collections - enter priority mode
      if ((ngcoll > 5) && (wm->GetNworking() <= 1)) {
        ngcoll = 0;
        for (Int_t slot = 0; slot < propagator->fNevents; slot++)
          if (propagator->fEvents[slot]->Prioritize())
            propagator->fPriorityEvents++;
        while ((!sch->GarbageCollect(td, true)) && (feederQ->size_async() == 0))
          ;
      }
      wm->FeederQueue()->wait_and_pop(basket);
    }
    // Check exit condition: null basket in the queue
    if (!basket)
      break;
    waiting[tid] = 0;
    if (td->NeedsToClean())
      td->CleanBaskets(0);
    else {
      // Check if there are at least 2 free baskets for this thread
      if (td->fPool.size() < 10) {
        if (basket->GetBasketMgr())
          basket->GetBasketMgr()->CreateEmptyBaskets(2, td);
        else
          td->fBmgr->CreateEmptyBaskets(10, td);
      }
    }
    ++counter;
    ntotransport = basket->GetNinput(); // all tracks to be transported
                                        //      ninput = ntotransport;
    GeantTrack_v &input = basket->GetInputTracks();
    GeantTrack_v &output = basket->GetOutputTracks();
    if (!ntotransport)
      goto finish; // input list empty
    //      Printf("======= BASKET %p with %d tracks counter=%d =======", basket, ntotransport,
    //      counter.load());
    //      basket->Print();
    //      Printf("==========================================");
    //      propagator->fTracksPerBasket[tid] = ntotransport;
    td->fVolume = 0;
    mat = 0;
    if (!basket->IsMixed()) {
      td->fVolume = basket->GetVolume();
#ifdef USE_VECGEOM_NAVIGATOR
      mat = ((TGeoMedium *)td->fVolume->getUserExtensionPtr())->GetMaterial();
#else
      mat = td->fVolume->GetMaterial();
#endif
      if (ntotransport < 257)
        nvect[ntotransport] += ntotransport;
    } else {
      nvect[1] += ntotransport;
    }

    // Select the discrete physics process for all particles in the basket
    if (propagator->fUsePhysics) {
      propagator->ProposeStep(ntotransport, input, td);
      // Apply msc for charged tracks
      propagator->ApplyMsc(ntotransport, input, td);
    }

    ncross = 0;
    generation = 0;

    while (ntotransport) {
      // Interrupt condition here. Work stealing could be also implemented here...
      generation++;
      //         Printf("====== WorkloadManager:");
      //         input.PrintTracks();
      // Propagate all remaining tracks
      if (basket->IsMixed())
        ncross += input.PropagateTracksScalar(output, td);
      else
        ncross += input.PropagateTracks(output, td);
      ntotransport = input.GetNtracks();
    }
    // All tracks are now in the output track vector. Possible statuses:
    // kCrossing - particles crossing boundaries
    // kPhysics - particles reaching the point where the discrete physics process
    //            will happen.
    // kExitingSetup - particles exiting the geometry
    // kKilled - particles that could not advance in geometry after several tries

    // Post-step actions by continuous processes for all particles. There are no
    // new generated particles at this point.
    if (propagator->fUsePhysics) {
      nphys = 0;
      nextra_at_rest = 0;
      // count phyics steps here
      Int_t nout = output.GetNtracks();
      for (auto itr = 0; itr < nout; ++itr)
        if (output.fStatusV[itr] == kPhysics)
          nphys++;
      if (nphys)
        propagator->fNphysSteps += nphys;

      propagator->Process()->Eloss(mat, output.GetNtracks(), output, nextra_at_rest, td);
      //         if (nextra_at_rest) Printf("Extra particles: %d", nextra_at_rest);
      // Now we may also have particles killed by energy threshold
      // Do post-step actions on remaining particles
      // to do: group particles per process

      if (propagator->fUsePhysics) {
        // Discrete processes only
        nphys = output.SortByLimitingDiscreteProcess(); // only those that kPhysics and not continous limit
        if (nphys) {
          // propagator->fNphysSteps += nphys;  dont do it here because dont count those killed in eloss
          // Do post step actions for particles suffering a given process.
          // Surviving particles are added to the output array

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
#if USE_VECPHYS == 1
          propagator->fVectorPhysicsProcess->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);
#endif
          // second: sample final states (based on the inf. regarding sampled
          //         target and type of interaction above), insert them into
          //         the track vector, update primary tracks;
          propagator->Process()->PostStepFinalStateSampling(mat, nphys, output, ntotnext, td);

          if (0 /*ntotnext*/) {
            printf("============= Basket: %s\n", basket->GetName());
            output.PrintTracks();
          }
        }
      }
    }
    if (gPropagator->fStdApplication)
      gPropagator->fStdApplication->StepManager(output.GetNtracks(), output, td);
    gPropagator->fApplication->StepManager(output.GetNtracks(), output, td);
    // Update geometry path for crossing tracks
    ntotnext = output.GetNtracks();

// Normalize directions (should be not needed)
//    for (auto itr=0; itr<ntotnext; ++itr)
//      output.Normalize(itr);

#ifdef BUG_HUNT
    for (auto itr = 0; itr < ntotnext; ++itr) {
      bool valid = output.CheckNavConsistency(itr);
      if (!valid) {
        valid = true;
      }
    }
    // First breakpoint to be set
    output.BreakOnStep(propagator->fDebugEvt, propagator->fDebugTrk, propagator->fDebugStp, prop->fDebugRep, "EndStep");
#endif
    for (auto itr = 0; itr < ntotnext; ++itr) {
      output.fNstepsV[itr]++;
      if (output.fStatusV[itr] == kBoundary)
        *output.fPathV[itr] = *output.fNextpathV[itr];
    }
  finish:
    //      basket->Clear();
    //      Printf("======= BASKET(tid=%d): in=%d out=%d =======", tid, ninput, basket->GetNoutput());
    //      ninjected =
    sch->AddTracks(basket, ntot, nnew, nkilled, td);
    //      Printf("thread %d: injected %d baskets", tid, ninjected);
    //      wm->TransportedQueue()->push(basket);
    //    sched_locker.StartOne();
    // Make sure the basket is not recycled before gets released by basketizer
    // This should not happen vey often, just when some threads are highly
    // demoted ant the basket makes it through the whole cycle before being fully released
    while (basket->fNused.load())
      ;
    basket->Recycle(td);
  }
  wm->DoneQueue()->push(0);
  delete prioritizer;
  Printf("=== Thread %d: exiting ===", tid);
  return 0;
}

//______________________________________________________________________________
void *WorkloadManager::TransportTracksCoprocessor(void *arg) {
  // Thread propagating all tracks from a basket.
  //      char slist[256];
  //      TString sslist;
  //   const Int_t max_idle = 1;
  //   Int_t indmin, indmax;
  static std::atomic<int> counter(0);
  Int_t ntotnext;
  // Int_t ncross;
  Int_t nbaskets;
  Int_t ntotransport;
  // Int_t nextra_at_rest = 0;
  Int_t generation = 0;
  //   Int_t ninjected = 0;
  Int_t nnew = 0;
  Int_t ntot = 0;
  Int_t nkilled = 0;
  Int_t ngcoll = 0;
  GeantBasket *basket = 0;

  Int_t tid = Instance()->ThreadId();
  Printf("=== Worker thread %d created for Coprocessor ===", tid);

  GeantPropagator *propagator = GeantPropagator::Instance();
  GeantTaskData *td = propagator->fThreadData[tid];
  td->fTid = tid;
  Int_t nworkers = propagator->fNthreads;
  WorkloadManager *wm = WorkloadManager::Instance();
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  GeantScheduler *sch = wm->GetScheduler();
  Int_t *nvect = sch->GetNvect();
  GeantBasketMgr *prioritizer = nullptr; // new GeantBasketMgr(sch, 0, 0);
  td->fBmgr = nullptr;                   // prioritizer;
  // prioritizer->SetThreshold(propagator->fNperBasket);
  // prioritizer->SetFeederQueue(wm->FeederQueue());
  // Start the feeder
  propagator->Feeder(td);
  // TGeoMaterial *mat = 0;
  Int_t *waiting = wm->GetWaiting();
//  condition_locker &sched_locker = wm->GetSchLocker();
// Int_t nprocesses = propagator->fNprocesses;
// Int_t ninput;
// Int_t noutput;
//   Bool_t useDebug = propagator->fUseDebug;
//   Printf("(%d) WORKER started", tid);
#ifdef USE_VECGEOM_NAVIGATOR
// Suppose I do not need a navigator, otherwise how it would have ever worked?
#else
  // Create navigator if none serving this thread.
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif
  waiting[tid] = 1;
  // Int_t iev[500], itrack[500];
  // TGeoBranchArray *crt[500], *nxt[500];

  TaskBroker *broker = reinterpret_cast<TaskBroker *>(arg);
  // broker->SetPrioritizer(prioritizer);
  while (1) {

    // Call the feeder if in priority mode
    if (prioritizer && // Disable this case until later.
        !prioritizer->HasTracks() && (propagator->GetNpriority() || wm->GetNworking() == 1)) {
      if (propagator->Feeder(td))
        ngcoll = 0;
      // Check exit condition
      if (propagator->TransportCompleted()) {
        for (Int_t i = 0; i < nworkers; i++)
          wm->FeederQueue()->push(0);
        wm->TransportedQueue()->push(0);
        wm->Stop();
        //         sched_locker.StartOne(); // signal the scheduler who has to exit
        //         gbc_locker.StartOne();
        break;
      }
    }

    // ::Info("GPU","Waiting (1) for next available stream.");
    TaskBroker::Stream stream = broker->GetNextStream();
    if (!stream)
      break;

    if (wm->FeederQueue()->empty()) {
      // There is no work to be done for now, let's just run what we have
      if (0 != broker->launchTask()) {
        // We ran something, let wait for the next free stream,
        continue;
      } else {
        // We had nothing to run at all .... we need to wait for
        // more data ....
      }
    }
    waiting[tid] = 1;
    nbaskets = feederQ->size_async();
    if (nbaskets > nworkers)
      ngcoll = 0;
    // If prioritizer has work, just do it
    if (prioritizer && prioritizer->HasTracks()) {
      basket = prioritizer->GetBasketForTransport(td);
      ngcoll = 0;
    } else {
      if (nbaskets < 1) {
        sch->GarbageCollect(td);
        ngcoll++;
      }
      // Too many garbage collections - enter priority mode
      if ((ngcoll > 5) && (wm->GetNworking() <= 1)) {
        ngcoll = 0;
        for (Int_t slot = 0; slot < propagator->fNevents; slot++)
          if (propagator->fEvents[slot]->Prioritize())
            propagator->fPriorityEvents++;
        while ((!sch->GarbageCollect(td, true)) && (feederQ->size_async() == 0))
          ;
      }
      wm->FeederQueue()->wait_and_pop(basket);
    }
    waiting[tid] = 0;
    if (td->NeedsToClean())
      td->CleanBaskets(0);
    else {
      // Check if there are at least 2 free baskets for this thread
      if (td->fPool.size() < 10) {
        if (basket->GetBasketMgr())
          basket->GetBasketMgr()->CreateEmptyBaskets(2, td);
        else
          td->fBmgr->CreateEmptyBaskets(10, td);
      }
    }
    // Check exit condition: null basket in the queue
    if (!basket)
      break;
    ++counter;

    if (!stream) {
      ::Info("GPU", "Waiting (2) for next available stream.");
      stream = broker->GetNextStream();
    }
    // lastToClear = false;
    if (!basket) {
      if (0 != broker->launchTask(/* wait= */ true)) {
        // We ran something, new basket might be available.
        continue;
      } else {
        break;
      }
    }

    ntotransport = basket->GetNinput(); // all tracks to be transported
    // ninput = ntotransport;
    GeantTrack_v &input = basket->GetInputTracks();
    GeantTrack_v &output = basket->GetOutputTracks();
    if (!ntotransport)
      goto finish; // input list empty
    //      Printf("======= BASKET %p with %d tracks counter=%d =======", basket, ntotransport,
    //      counter);
    //      basket->Print();
    //      Printf("==========================================");
    //      propagator->fTracksPerBasket[tid] = ntotransport;
    td->fVolume = 0;
    // mat = 0;
    if (!basket->IsMixed()) {
      td->fVolume = basket->GetVolume();
      // mat = td->fVolume->GetMaterial();
      if (ntotransport < 257)
        nvect[ntotransport] += ntotransport;
    } else {
      nvect[1] += ntotransport;
    }

    // Record tracks
    // ninput = ntotransport;
    if (basket->GetNoutput()) {
      Geant::Warning("TransportTracksCoprocessor", "Output Track_v not empty noutput=%d counter=%d",
                     basket->GetNoutput(), counter.load());
    }
    //      if (counter==1) input.PrintTracks();
    for (Int_t itr = 0; itr < ntotransport; itr++) {
      // iev[itr] = input.fEventV[itr];
      // itrack[itr] = input.fParticleV[itr];
      // crt[itr] = input.fPathV[itr];
      // nxt[itr] = input.fNextpathV[itr];
      if (TMath::IsNaN(input.fXdirV[itr])) {
        Printf("Error: track %d has NaN", itr);
      }
    }
    // Select the discrete physics process for all particles in the basket
    // if (propagator->fUsePhysics) {
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
      // NOTE: need to deal with propagator->fUsePhysics
      broker->runTask(tid, *basket); // ntotransport, basket_sch->GetNumber(), gPropagator->fTracks, particles);
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

    // Check
    if (basket->GetNinput()) {
      Geant::Warning("TransportTracksCoprocessor", "Input Track_v not empty: ninput=%d noutput=%d counter=%d",
                     basket->GetNinput(), basket->GetNoutput(), counter.load());
    }
    {
      auto noutput = basket->GetNoutput();
      for (Int_t itr = 0; itr < noutput; itr++) {
        if (TMath::IsNaN(output.fXdirV[itr])) {
          Geant::Error("TransportTracksCoprocessor", "Track %d has NaN", itr);
        }
      }
    }
    if (gPropagator->fStdApplication)
      gPropagator->fStdApplication->StepManager(output.GetNtracks(), output, td);
    gPropagator->fApplication->StepManager(output.GetNtracks(), output, td);
    // Update geometry path for crossing tracks
    ntotnext = output.GetNtracks();

#ifdef BUG_HUNT
    for (auto itr = 0; itr < ntotnext; ++itr) {
      bool valid = output.CheckNavConsistency(itr);
      if (!valid) {
        valid = true;
      }
    }
    // First breakpoint to be set
    output.BreakOnStep(propagator->fDebugEvt, propagator->fDebugTrk, propagator->fDebugStp, prop->fDebugRep, "EndStep");
#endif
    for (auto itr = 0; itr < ntotnext; ++itr) {
      output.fNstepsV[itr]++;
      if (output.fStatusV[itr] == kBoundary)
        *output.fPathV[itr] = *output.fNextpathV[itr];
    }
  finish:
    //      basket->Clear();
    //      Printf("======= BASKET(tid=%d): in=%d out=%d =======", tid, ninput, basket->GetNoutput());
    /* Int_t ninjected = */ sch->AddTracks(basket, ntot, nnew, nkilled, nullptr /* prioritizer */);
    //      Printf("thread %d: injected %d baskets", tid, ninjected);
    // wm->TransportedQueue()->push(basket);
    (void)ntot;
    (void)nnew;
    (void)nkilled;
    //    sched_locker.StartOne();
    basket->Recycle(td);
  }
  // delete prioritizer;
  wm->DoneQueue()->push(0);
  Printf("=== Coprocessor Thread %d: exiting === Processed %ld", tid, broker->GetTotalWork());
  return 0;
}

//______________________________________________________________________________
void *WorkloadManager::GarbageCollectorThread(void *) {
  // This threads can be triggered to do garbage collection of unused baskets
  static Double_t rsslast = 0;
  Double_t rss;
  ProcInfo_t procInfo;
  const Double_t MByte = 1024.;
  const double_t thr_increase = 1.05;
  GeantPropagator *propagator = GeantPropagator::Instance();
  WorkloadManager *wm = WorkloadManager::Instance();
  Int_t nthreads = propagator->fNthreads;
  Double_t threshold = propagator->fMaxRes;
  if (threshold == 0)
    return 0;
  //  condition_locker &gbc_locker = wm->GetGbcLocker();
  while (1) {
    //    gbc_locker.Wait();
    if (wm->IsStopped())
      break;
    // todo: cleaning with thread based recycle queues
    gSystem->GetProcInfo(&procInfo);
    rss = procInfo.fMemResident / MByte;
    //    Printf("### Checking mem: %g MBytes", rss);
    if (rss > threshold && (rss / rsslast > thr_increase)) {
      rsslast = rss;
      for (Int_t tid = 0; tid < nthreads; tid++) {
        GeantTaskData *td = propagator->fThreadData[tid];
        td->SetToClean(true);
      }
    }
    gSystem->Sleep(1000); // millisec
  }
  return 0;
}

//______________________________________________________________________________
Int_t WorkloadManager::GetMonFeatures() const {
  // Get the number of monitored features
  return (fMonQueue + fMonMemory + fMonBasketsPerVol + fMonVectors + fMonConcurrency + fMonTracksPerEvent + fMonTracks);
}

//______________________________________________________________________________
bool WorkloadManager::IsMonitored(EGeantMonitoringType feature) const {
  // Check if a given feature is monitored
  switch (feature) {
  case kMonQueue:
    return fMonQueue;
  case kMonMemory:
    return fMonMemory;
  case kMonBasketsPerVol:
    return fMonBasketsPerVol;
  case kMonVectors:
    return fMonVectors;
  case kMonConcurrency:
    return fMonConcurrency;
  case kMonTracksPerEvent:
    return fMonTracksPerEvent;
  case kMonTracks:
    return fMonTracks;
  }
  return false;
}

//______________________________________________________________________________
void WorkloadManager::SetMonitored(EGeantMonitoringType feature, bool flag) {
  // Enable/disable monitoring for a feature
  int value = (int)flag;
  switch (feature) {
  case kMonQueue:
    fMonQueue = value;
    break;
  case kMonMemory:
    fMonMemory = value;
    break;
  case kMonBasketsPerVol:
    fMonBasketsPerVol = value;
    break;
  case kMonVectors:
    fMonVectors = value;
    break;
  case kMonConcurrency:
    fMonConcurrency = value;
    break;
  case kMonTracksPerEvent:
    fMonTracksPerEvent = value;
    break;
  case kMonTracks:
    fMonTracks = value;
  }
}

//______________________________________________________________________________
void *WorkloadManager::MonitoringThread(void *) {
  // Thread providing basic monitoring for the scheduler.
  const Double_t MByte = 1024.;
  Printf("Started monitoring thread...");
  GeantPropagator *propagator = GeantPropagator::Instance();
  WorkloadManager *wm = WorkloadManager::Instance();
  Int_t nmon = wm->GetMonFeatures();
  if (!nmon)
    return 0;
  Int_t dmon = 0.5 * nmon + nmon % 2;
  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  Int_t ntotransport;
  ProcInfo_t procInfo;
  Double_t rss;
  Double_t nmem[101] = {0};
  GeantScheduler *sch = wm->GetScheduler();
  Int_t nvol = sch->GetNvolumes();
  Int_t *nvect = sch->GetNvect();
  Int_t nthreads = wm->GetNthreads();
  Int_t *nworking = new Int_t[nthreads + 1];
  memset(nworking, 0, (nthreads + 1) * sizeof(Int_t));
  Int_t *waiting = wm->GetWaiting();
  GeantBasketMgr **bmgr = sch->GetBasketManagers();

  TCanvas *cmon = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("cscheduler");
  if (nmon == 1)
    cmon->Divide(1, 1);
  else
    cmon->Divide(2, dmon);
  TH1I *hqueue = 0;
  Int_t nqueue[101] = {0};
  Int_t ipad = 0;
  if (wm->IsMonitored(WorkloadManager::kMonQueue)) {
    hqueue = new TH1I("hqueue", "Work queue load", 100, 0, 100);
    hqueue->SetFillColor(kRed);
    hqueue->SetLineColor(0);
    hqueue->SetStats(false);
    cmon->cd(++ipad);
    hqueue->Draw();
  }
  TH1F *hmem = 0;
  if (wm->IsMonitored(WorkloadManager::kMonMemory)) {
    hmem = new TH1F("hmem", "Resident memory [MB]", 100, 0, 100);
    // hmem->SetFillColor(kMagenta);
    hmem->SetLineColor(kMagenta);
    hmem->SetStats(false);
    cmon->cd(++ipad);
    hmem->Draw();
  }
  TH1I *hbaskets = 0;
  TH1I *hbused = 0;
  if (wm->IsMonitored(WorkloadManager::kMonBasketsPerVol)) {
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
  if (wm->IsMonitored(WorkloadManager::kMonVectors)) {
    hvectors = new TH1I("hvectors", "Tracks in vectors of given size", 257, 0, 257);
    hvectors->SetFillColor(kBlue);
    hvectors->SetLineColor(0);
    hvectors->SetStats(false);
    cmon->cd(++ipad)->SetLogy();
    hvectors->Draw();
  }
  TH1F *hconcurrency = 0;
  TH1F *hconcavg = 0;
  if (wm->IsMonitored(WorkloadManager::kMonConcurrency)) {
    hconcurrency = new TH1F("hconcurrency", "Concurrency plot", nthreads + 1, 0, nthreads + 1);
    hconcurrency->GetYaxis()->SetRangeUser(0, 1);
    hconcurrency->GetXaxis()->SetNdivisions(nthreads + 1, true);
    hconcurrency->GetYaxis()->SetNdivisions(10, true);
    hconcurrency->SetFillColor(kGreen);
    hconcurrency->SetLineColor(0);
    hconcurrency->SetStats(false);
    hconcavg = new TH1F("hconcavg", "Concurrency plot", nthreads + 1, 0, nthreads + 1);
    hconcavg->SetFillColor(kRed);
    hconcavg->SetFillStyle(3001);
    hconcavg->SetLineColor(0);
    hconcavg->SetStats(false);
    cmon->cd(++ipad);
    hconcurrency->Draw();
    hconcavg->Draw("SAME");
  }
  TH1I *htracksmax = 0;
  TH1I *htracks = 0;
  Int_t nbuffered = propagator->fNevents;
  if (wm->IsMonitored(WorkloadManager::kMonTracksPerEvent)) {
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
  Int_t ntrackstot[101] = {0};
  if (wm->IsMonitored(WorkloadManager::kMonTracks)) {
    htrackstot = new TH1I("htrackstot", "Total number of tracks alive", 100, 0, 100);
    htrackstot->SetFillColor(kRed);
    htrackstot->SetLineColor(0);
    htrackstot->SetStats(false);
    cmon->cd(++ipad);
    htrackstot->Draw();
  }
  cmon->Update();
  Double_t stamp = 0.;
  Int_t i, j, bin;
  Int_t nmaxtot;
  while (1) { // exit condition here
    i = Int_t(stamp);
    ipad = 0;
    gSystem->Sleep(50); // millisec
    // Fill histograms
    if (stamp > 100) {
      if (hqueue) {
        ntotransport = feederQ->size_async();
        memmove(nqueue, &nqueue[1], 99 * sizeof(Int_t));
        nqueue[99] = ntotransport;
        hqueue->GetXaxis()->Set(100, stamp - 100, stamp);
        for (j = 0; j < 100; j++)
          hqueue->SetBinContent(j + 1, nqueue[j]);
      }
      if (hmem) {
        gSystem->GetProcInfo(&procInfo);
        rss = procInfo.fMemResident / MByte;
        memmove(nmem, &nmem[1], 99 * sizeof(Double_t));
        nmem[99] = rss;
        hmem->GetXaxis()->Set(100, stamp - 100, stamp);
        for (j = 0; j < 100; j++)
          hmem->SetBinContent(j + 1, nmem[j]);
      }
      if (htrackstot) {
        // Count tracks for all event slots
        Int_t ntr = 0;
        for (Int_t slot = 0; slot < propagator->fNevents; slot++)
          ntr += propagator->fEvents[slot]->GetNinflight();
        memmove(ntrackstot, &ntrackstot[1], 99 * sizeof(Int_t));
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
        nmem[i] = rss;
        hmem->SetBinContent(i + 1, nmem[i]);
      }
      if (htrackstot) {
        // Count tracks for all event slots
        Int_t ntr = 0;
        for (Int_t slot = 0; slot < propagator->fNevents; slot++)
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
        GeantTaskData *td = propagator->fThreadData[j];
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
    if (hconcurrency) {
      for (j = 0; j < nthreads + 1; j++) {
        nworking[j] += 1 - waiting[j];
        hconcurrency->SetBinContent(j + 1, 1 - waiting[j]);
        hconcavg->SetBinContent(j + 1, nworking[j] / (stamp + 1));
      }
    }
    if (htracksmax) {
      nmaxtot = 1;
      for (j = 0; j < nbuffered; j++) {
        GeantEvent *evt = propagator->fEvents[j];
        Int_t nmax = evt->GetNmax();
        nmaxtot = TMath::Max(nmax, nmaxtot);
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
    if (hconcurrency) {
      cmon->cd(++ipad);
      hconcurrency->Draw();
      hconcavg->Draw("SAME");
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
  Printf("=== Percent of tracks transported in single track mode: %g%%", 100. * nvect[1] / sumvect);
  // Sleep a bit to let the graphics finish
  gSystem->Sleep(100); // millisec

  return 0;
}
