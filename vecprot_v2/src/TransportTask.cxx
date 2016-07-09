#include "Geant/Error.h"

#include "TransportTask.h"
#include "ThreadData.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "GeantEvent.h"

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
#include "GeantTrack.h"
#include "GeantBasket.h"
#include "GeantOutput.h"
#include "GeantTaskData.h"
#include "PhysicsProcess.h"
#include "GeantScheduler.h"
#include "GeantEvent.h"
#include "GeantVApplication.h"
#if USE_VECGEOM_NAVIGATOR
#include "base/TLS.h"
#include "management/GeoManager.h"
#include "materials/Medium.h"
#else
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#endif
#include "TaskBroker.h"

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

TransportTask::TransportTask () { }

TransportTask::~TransportTask () { }

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

tbb::task* TransportTask::execute ()
{
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
  //   int ninjected = 0;
  int nnew = 0;
  int ntot = 0;
  int nkilled = 0;
  int nphys = 0;
  int nout = 0;
  int ngcoll = 0;
  GeantBasket *basket = 0;
  GeantPropagator *propagator = GeantPropagator::Instance();
  ThreadData *threadData = ThreadData::Instance(propagator->fNthreads);
  WorkloadManager *wm = WorkloadManager::Instance();
  //if (wm->IsStopped()){
  //  tbb::task::cancel_group_execution();
  //  return NULL;
  //}
  int tid = wm->ThreadId();
  Geant::Print("","=== Worker thread %d created ===", tid);
  if(tid>=propagator->fNthreads){
    Error("TransportTask", "Task was put in a non-working thread ");
    return NULL;
  }
  Geant::GeantTaskData *td = propagator->fThreadData[tid];
  td->fTid = tid;
  int nworkers = propagator->fNthreads;

  Geant::priority_queue<GeantBasket *> *feederQ = wm->FeederQueue();
  GeantScheduler *sch = wm->GetScheduler();
  int *nvect = sch->GetNvect();
  GeantBasketMgr *prioritizer = threadData->fPrioritizers[tid];
  td->fBmgr = prioritizer;
  prioritizer->SetThreshold(propagator->fNperBasket);
  prioritizer->SetFeederQueue(feederQ);

  // IO handling

  bool concurrentWrite = GeantPropagator::Instance()->fConcurrentWrite && GeantPropagator::Instance()->fFillTree;
  int treeSizeWriteThreshold = GeantPropagator::Instance()->fTreeSizeWriteThreshold;

  GeantFactory<MyHit> *myhitFactory = threadData->fMyhitFactories[tid];

  TThread t;
  TThreadMergingFile* file = threadData->fFiles[tid];
  TTree *tree = threadData->fTrees[tid];
  GeantBlock<MyHit>* data = threadData->fData[tid];

  /*
  if (concurrentWrite)
    {
      file = new TThreadMergingFile("hits_output.root", wm->IOQueue(), "RECREATE");
      tree = new TTree("Tree","Simulation output");

      tree->Branch("hitblockoutput", "GeantBlock<MyHit>", &data);

      // set factory to use thread-local queues
      myhitFactory->queue_per_thread = true;
    }*/

  // Start the feeder
  //propagator->Feeder(td);
  /*int returning;
  FeederTask & feederTask = *new(tbb::task::allocate_root()) FeederTask(td, &returning);
  tbb::task::spawn_root_and_wait(feederTask); */

  Material_t *mat = 0;
  int *waiting = wm->GetWaiting();
//  condition_locker &sched_locker = wm->GetSchLocker();
//  condition_locker &gbc_locker = wm->GetGbcLocker();
// The last worker wakes the scheduler
//  if (tid == nworkers-1) sched_locker.StartOne();
//   int nprocesses = propagator->fNprocesses;
//   int ninput, noutput;
//   bool useDebug = propagator->fUseDebug;
//   Geant::Print("","(%d) WORKER started", tid);
// Create navigator if none serving this thread.
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
    auto feedres = wm->CheckFeederAndExit(*prioritizer, *propagator, *td);
    if (feedres == WorkloadManager::FeederResult::kFeederWork) {
       ngcoll = 0;
    } else if (feedres == WorkloadManager::FeederResult::kStopProcessing) {
       break;
    }

    // Collect info about the queue
    waiting[tid] = 1;
    nbaskets = feederQ->size_async();
    if (nbaskets > nworkers)
      ngcoll = 0;

    // Fire garbage collection if starving
    if ((nbaskets < 1) && (!propagator->IsFeeding())) {
      sch->GarbageCollect(td);
     ngcoll++;
    }
    // Too many garbage collections - enter priority mode
    if ((ngcoll > 5) && (wm->GetNworking() <= 1)) {
      ngcoll = 0;
      for (int slot = 0; slot < propagator->fNevents; slot++)
        if (propagator->fEvents[slot]->Prioritize())
          propagator->fPriorityEvents++;
      while ((!sch->GarbageCollect(td, true)) &&
             (feederQ->size_async() == 0) &&
             (!basket) &&
             (!prioritizer->HasTracks()))
        ;
    }
    // Check if the current basket is reused or we need a new one
    if (!basket) {
      // If prioritizer has work, just do it
      if (prioritizer->HasTracks()) {
        basket = prioritizer->GetBasketForTransport(td);
        ngcoll = 0;
      } else {
        // Take next basket from queue
        wm->FeederQueue()->wait_and_pop(basket);
        // If basket from queue is null, exit
        if (!basket)
          break;
      }
    }
    // Start transporting the basket
    waiting[tid] = 0;
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
      mat = ((Medium_t *)td->fVolume->GetTrackingMediumPtr())->GetMaterial();
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
      // Use fNsteps track data to detect geometry anomalies
      for (auto itr=0; itr<ntotransport; ++itr) {
        input.fNstepsV[itr]++;
        if ((input.fStatusV[itr] != kKilled) &&
            (input.fNstepsV[itr] > gPropagator->fNstepsKillThr) &&
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
    if (propagator->fUsePhysics) {
      nphys = 0;
      nextra_at_rest = 0;
      // count phyics steps here
      for (auto itr = 0; itr < nout; ++itr)
        if (output.fStatusV[itr] == kPhysics)
          nphys++;

      propagator->Process()->Eloss(mat, output.GetNtracks(), output, nextra_at_rest, td);
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

        if (0 /*ntotnext*/) {
          Geant::Print("","============= Basket: %s\n", basket->GetName());
          output.PrintTracks();
        }
      }
    }
    if (gPropagator->fStdApplication)
      gPropagator->fStdApplication->StepManager(output.GetNtracks(), output, td);
    gPropagator->fApplication->StepManager(output.GetNtracks(), output, td);

    // WP
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
    output.BreakOnStep(propagator->fDebugEvt, propagator->fDebugTrk, propagator->fDebugStp, propagator->fDebugRep, "EndStep");
#endif
    for (auto itr = 0; itr < ntotnext; ++itr) {
      //output.fNstepsV[itr]++;
      if (output.fStatusV[itr] == kBoundary)
        *output.fPathV[itr] = *output.fNextpathV[itr];
    }
  finish:
    // Check if there are enough transported tracks staying in the same volume
    // to be reused without re-basketizing
    int nreusable = sch->ReusableTracks(output);
    bool reusable = (basket->IsMixed())? false : (nreusable>=gPropagator->fNminReuse);
//    reusable = false;

    if (reusable)
      sch->CopyReusableTracks(output, input, basket->GetThreshold());
    // Remaining tracks need to be re-basketized
    sch->AddTracks(output, ntot, nnew, nkilled, td);
    // Make sure the basket is not recycled before gets released by basketizer
    // This should not happen vey often, just when some threads are highly
    // demoted ant the basket makes it through the whole cycle before being fully released
    if (!reusable) {
      // Recycle basket, otherwise keep it for the next iteration
      while (basket->fNused.load())
        ;
      basket->Recycle(td);
      basket = 0; // signal reusability by non-null pointer
    }
    // Update boundary crossing counter
    td->fNcross += ncross;
  }

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

  Geant::Print("","=== Thread %d: exiting ===", tid);

  if (wm->IsStopped()) wm->MergingServer()->Finish();

  if (concurrentWrite) {
    delete file;
  }

  return NULL;
}
