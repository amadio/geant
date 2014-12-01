#include "WorkloadManager.h"

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
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "GeantTrack.h"
#include "GeantBasket.h"
#include "GeantOutput.h"
#include "GeantThreadData.h"
#include "PhysicsProcess.h"
#include "GeantScheduler.h"
#include "GeantEvent.h"
#include "GeantVApplication.h"
#if USE_VECGEOM_NAVIGATOR == 1
 #include "management/GeoManager.h"
#endif
#include "TaskBroker.h"

ClassImp(WorkloadManager)

WorkloadManager *WorkloadManager::fgInstance = 0;

//______________________________________________________________________________
WorkloadManager::WorkloadManager(Int_t nthreads)
                :TObject(),
                 fNthreads(nthreads),
                 fNbaskets(0),
                 fBasketGeneration(0),
                 fNbasketgen(0),
                 fNidle(nthreads),
                 fNminThreshold(10),
                 fNqueued(0),
                 fBtogo(0),
                 fStarted(kFALSE),
                 fStopped(kFALSE),
                 fFeederQ(0),
                 fTransportedQ(0),
                 fDoneQ(0),
//                 fNavStates(0),
                 fListThreads(0),
                 fFlushed(kFALSE),
                 fFilling(kFALSE),
                 fScheduler(0),
                 fBroker(0),
                 fWaiting(0),
                 fMutexSch(new std::mutex),
                 fCondSch(new std::condition_variable),
                 fWorkDone(false)
{
// Private constructor.
   fFeederQ = new Geant::priority_queue<GeantBasket*>(1<<16);
   fTransportedQ = new Geant::priority_queue<GeantBasket*>(1<<16);
   fDoneQ = new Geant::priority_queue<GeantBasket*>(1<<10);
   fgInstance = this;
   fScheduler = new GeantScheduler();
   fWaiting = new Int_t[fNthreads+1];
   memset(fWaiting, 0, (fNthreads+1)*sizeof(Int_t));   
}

//______________________________________________________________________________
WorkloadManager::~WorkloadManager()
{
// Destructor.
   delete fFeederQ;
   delete fTransportedQ;
   delete fDoneQ;
   delete fScheduler;
   delete [] fWaiting;
//   delete fNavStates;
   fgInstance = 0;
}

//______________________________________________________________________________
void WorkloadManager::CreateBaskets()
{
// Create the array of baskets
   VolumePath_t *blueprint = 0;
   Int_t maxdepth = GeantPropagator::Instance()->fMaxDepth;
   Printf("Max depth: %d", maxdepth);
   blueprint = VolumePath_t::MakeInstance(maxdepth);
//   fNavStates = new GeantObjectPool<VolumePath_t>(1000*fNthreads, blueprint);
//   fNavStates = new rr_pool<VolumePath_t>(16*fNthreads, 1000, blueprint);
   fScheduler->CreateBaskets();
   VolumePath_t::ReleaseInstance(blueprint);
}
   
//______________________________________________________________________________
WorkloadManager *WorkloadManager::Instance(Int_t nthreads)
{
// Return singleton instance.
   if (fgInstance) return fgInstance;
   if (!nthreads) {
      ::Error("WorkloadManager::Instance", "No instance yet so you should provide number of threads.");
      return 0;
   }
   return new WorkloadManager(nthreads);
}

//______________________________________________________________________________
void WorkloadManager::Print(Option_t *) const
{
//
}   

//______________________________________________________________________________
void WorkloadManager::SetTaskBroker(TaskBroker *broker) { 
  // Register the broker if it is valid.
  if (broker && broker->IsValid()) fBroker = broker; 
}


//______________________________________________________________________________
void WorkloadManager::StartThreads()
{
// Start the threads
   fStarted = kTRUE;
   if (fListThreads) return;
   fListThreads = new TList();
   fListThreads->SetOwner();
   Int_t ith = 0;
   if (fBroker) {
     if (fBroker->GetNstream() > fNthreads) {
        Fatal("StartThreads","The task broker is using too many threads (%d out of %d)",
              fBroker->GetNstream(), fNthreads);
     }
     Printf("Running with a coprocessor broker.");
     TThread *t = new TThread(WorkloadManager::TransportTracksCoprocessor,fBroker);
     fListThreads->Add(t);
     t->Run();
     ith += fBroker->GetNstream() + 1;
   }
   for (; ith<fNthreads; ith++) {
      TThread *t = new TThread(WorkloadManager::TransportTracks);
      fListThreads->Add(t);
      t->Run();
   }
//   gSystem->Sleep(1000);
   TThread *t = new TThread(WorkloadManager::MainScheduler);
   fListThreads->Add(t);
   t->Run();
   if (GeantPropagator::Instance()->fUseMonitoring) {
      t = new TThread(WorkloadManager::MonitoringThread);
      fListThreads->Add(t);
      t->Run();
   }   
}   

//______________________________________________________________________________
void WorkloadManager::JoinThreads()
{
// 
   int tojoin = fNthreads;
   if (fBroker) tojoin -= fBroker->GetNstream();
   for (Int_t ith=0; ith<tojoin; ith++) fFeederQ->push(0);
   for (Int_t ith=0; ith<tojoin; ith++) ((TThread*)fListThreads->At(ith))->Join();
   // Join garbage collector
   ((TThread*)fListThreads->At(tojoin))->Join();
   // Join monitoring thread
   if (GeantPropagator::Instance()->fUseMonitoring) {
      ((TThread*)fListThreads->At(tojoin+1))->Join();
   }   
}
   
//______________________________________________________________________________
void WorkloadManager::WaitWorkers()
{
// Waiting point for the main thread until work gets done.
   fFilling = kFALSE;
   Int_t ntowait = fNthreads;
   if (fBroker) ntowait -= fBroker->GetNstream();
   GeantBasket *signal;
   while(ntowait) {
      fDoneQ->wait_and_pop(signal);
      ntowait--;
      Printf("=== %d workers finished", fNthreads+1-ntowait);
   }
//   fBasketGeneration++;
}

//______________________________________________________________________________
void *WorkloadManager::MainScheduler(void *)
{
// Garbage collector thread, called by a single thread.
   GeantPropagator *propagator = GeantPropagator::Instance();
   
   Int_t nworkers = propagator->fNthreads;
//   Int_t naverage = propagator->fNaverage;
//   Int_t maxperevent = propagator->fMaxPerEvent;
   Int_t dumped_event = -1;
   Int_t first_not_transported = 0;
   Int_t nbuffered = propagator->fNevents;
   Int_t last_event = nbuffered;
   Int_t max_events = propagator->fNtotal;
   TBits finished(max_events);
   Bool_t prioritize = kFALSE;
   Bool_t countdown = kFALSE;
   WorkloadManager *wm = WorkloadManager::Instance();
   if (wm->fBroker) nworkers -= wm->fBroker->GetNstream();
   Geant::priority_queue<GeantBasket*> *feederQ = wm->FeederQueue();
   //priority_queue<GeantBasket*> *transportedQ = wm->TransportedQueue();
   GeantScheduler *sch = wm->GetScheduler();
   // Number of baskets in the queue to transport
   Int_t ntotransport = 0;
   Int_t *waiting = wm->GetWaiting();
//   Int_t ntracksperbasket = propagator->fNperBasket;
   // Feeder threshold
   // Model parameters here
   Int_t min_feeder = TMath::Max(20,3*nworkers); // setter/getter here ?
   Int_t max_feeder = 300; // ???
   // Number of tracks in the current basket
//   Int_t ntracksb = 0;
   // Number of priority baskets
   Int_t npriority;
   //   Int_t lastphase = -1;
   Int_t crtphase  = 0;
   
   Int_t niter = -1;
//   size_t npop = 0;
   Int_t ninjected = 0;
//   Int_t nnew = 0, ntot=0, nkilled=0;
//   GeantBasket *output;
//   GeantBasket **carray = new GeantBasket*[500];
   waiting[nworkers] = 1;
   std::mutex *mutex_sch = wm->GetMutexSch();
   std::condition_variable *cond_sch = wm->GetCondSch();
   while (1) {
      std::unique_lock<std::mutex> lk(*mutex_sch);
      while (!wm->IsWorkDone()) cond_sch->wait(lk);
      wm->SetWorkDone(false);
      lk.unlock();
//      transportedQ->wait_and_pop_max(500,npop,carray);
      waiting[nworkers] = 0;
      // Monitor the queues while there are tracks to transport
      niter++;
//      if (niter==1000) exit(0);
      ninjected = 0;
//      nnew = 0;
//      ntot = 0;
//      nkilled = 0;
//      for (size_t iout=0; iout<npop; iout++) {
//         output = carray[iout];
//         ninjected += sch->AddTracks(output, ntot, nnew, nkilled);
         // Recycle basket
//         output->Recycle();
//  }
      // If there were events to be dumped, check their status here
      ntotransport = feederQ->size_async();
//      Printf("#%d: feeder=%p Processed %d baskets (%d tracks, %d new, %d killed)-> injected %d. QS=%d", niter, feederQ, npop, ntot, nnew, nkilled, ninjected, ntotransport);
#ifdef __STAT_DEBUG
           sch->GetPendingStat().Print();
           sch->GetQueuedStat().Print();
           sch->GetTransportStat().Print();
#endif           
      
      // Check and mark finished events
      for (Int_t ievt=0; ievt<nbuffered; ievt++) {
         GeantEvent *evt = propagator->fEvents[ievt];
         if (prioritize && ievt==dumped_event) evt->Print();
         if (finished.TestBitNumber(evt->GetEvent())) continue;
         if (evt->Transported()) {
            // Digitizer (delete for now)
            Int_t ntracks = propagator->fNtracks[ievt];
            Printf("= digitizing event %d with %d tracks", evt->GetEvent(), ntracks);
//            propagator->fApplication->Digitize(evt->GetEvent());
//            for (Int_t itrack=0; itrack--------------------------------------------------<ntracks; itrack++) {
//               delete propagator->fTracks[maxperevent*ievt+itrack];
//               propagator->fTracks[maxperevent*ievt+itrack] = 0;
//            }
            finished.SetBitNumber(evt->GetEvent());
            if (last_event<max_events) {
               Printf("=> Importing event %d", last_event);
               ninjected += propagator->ImportTracks(1,propagator->fNaverage,last_event,ievt);
               last_event++;
            }
         }
      }
      // Exit condition
      if (finished.FirstNullBit() >= (UInt_t)max_events) break;
      
      // In case some events were transported with priority, check if they finished
      if (prioritize) {
         first_not_transported = finished.FirstNullBit();
         if (first_not_transported > dumped_event) {
            // Priority events digitized, exit prioritized regime
            Printf("= stopped prioritizing");
            prioritize = kFALSE;
            sch->SetPriorityRange(-1, -1);
         } else {
            // Flush priority baskets if the countdown is zero
//            npriority = sch->GetNpriority();
            if (countdown) {
               if (feederQ->get_countdown()==0)  {
                  countdown = kFALSE;
                  feederQ->reset_countdown();
                  npriority = sch->FlushPriorityBaskets();
                  ninjected += npriority;
//                  Printf("Flushed %d priority baskets, resetting countdown", npriority);
               } else {
//                  Printf("Countdown is %d", feederQ->get_countdown());
               }
            } else {
               npriority = sch->FlushPriorityBaskets();
               ninjected = npriority;
//               Printf("Flushed %d priority baskets", npriority);
#ifdef __STAT_DEBUG
           Printf("After FlushPriorityBaskets:");
           sch->GetPendingStat().Print();
           sch->GetQueuedStat().Print();
           sch->GetTransportStat().Print();
#endif           
            }
         }
      }
      
      ntotransport = feederQ->size_async();
      if (ntotransport<min_feeder || ntotransport>max_feeder) sch->AdjustBasketSize();
      if (ntotransport<min_feeder) {
     // Transport queue below the threshold
         if (crtphase<1) crtphase = 1;
        // Set the countdown to the number of remaining objects
        
        // In case no new events were injected and we are not in a the priority regime
        // and below lowest watermark in this iteration, make a garbage collection
//        if (!nnew && !prioritize && ntotransport<nworkers) {

//           ninjected += sch->GarbageCollect();
//           Printf("Garbage collection injected %d baskets", ninjected);

//        }  
        if (!prioritize && last_event<max_events) {
           // Start prioritized regime
           dumped_event = finished.FirstNullBit();
           sch->SetPriorityRange(dumped_event, dumped_event+4);
#ifdef __STAT_DEBUG
           Printf("Before CollectPrioritizedTracks:");
           sch->GetPendingStat().Print();
           sch->GetQueuedStat().Print();
           sch->GetTransportStat().Print();
#endif           
           ninjected += sch->CollectPrioritizedTracks();
#ifdef __STAT_DEBUG
           Printf("After CollectPrioritizedTracks:");
           sch->GetPendingStat().Print();
           sch->GetQueuedStat().Print();
           sch->GetTransportStat().Print();
#endif           
           prioritize = kTRUE;
           countdown = kTRUE;
           ntotransport = feederQ->size_async();
           feederQ->set_countdown(ntotransport);
           Printf("====== Prioritizing events %d to %d, countdown=%d", dumped_event,dumped_event+4, ntotransport);
           waiting[nworkers] = 1;
           continue;
        }
      }
      ntotransport = feederQ->size_async();
      if (ntotransport==0) {
         Printf("Garbage collection");
         sch->GarbageCollect();
         if (countdown) feederQ->set_countdown(0);
      }
      waiting[nworkers] = 1;
   }
   // Stop workers
   //gROOT->SetBit(TObject::kInvalidObject, kTRUE);
   for (Int_t i=0; i<nworkers; i++) wm->FeederQueue()->push(0);
   wm->TransportedQueue()->push(0);
   wm->Stop();
   //gROOT->SetBit(TObject::kInvalidObject, kFALSE);
   propagator->fApplication->Digitize(0);
   Printf("=== Scheduler: stopping threads and exiting === niter =%d\n", niter);
   return 0;
}        

//______________________________________________________________________________
void *WorkloadManager::TransportTracks(void *)
{
// Thread propagating all tracks from a basket.
//      char slist[256];
//      TString sslist;
//   const Int_t max_idle = 1;
//   Int_t indmin, indmax;
   static std::atomic<int> counter(0);
   Int_t ntotnext, ncross;
   Int_t ntotransport;
   Int_t nextra_at_rest = 0;
   Int_t generation = 0;
//   Int_t ninjected = 0;
   Int_t nnew = 0;
   Int_t ntot = 0;
   Int_t nkilled = 0;
   GeantBasket *basket = 0;
   Int_t tid = TGeoManager::ThreadId();
   GeantPropagator *propagator = GeantPropagator::Instance();
   GeantThreadData *td = propagator->fThreadData[tid];
   WorkloadManager *wm = WorkloadManager::Instance();
   GeantScheduler *sch = wm->GetScheduler();
   Int_t *waiting = wm->GetWaiting();
//   std::mutex *mutex_sch = wm->GetMutexSch();
   std::condition_variable *cond_sch = wm->GetCondSch();
//   Int_t nprocesses = propagator->fNprocesses;
//   Int_t ninput, noutput;
//   Bool_t useDebug = propagator->fUseDebug;
//   Printf("(%d) WORKER started", tid);
   // Create navigator if none serving this thread.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) nav = gGeoManager->AddNavigator();
//   Int_t iev[500], itrack[500];
   // TGeoBranchArray *crt[500], *nxt[500];
   while (1) {
      waiting[tid] = 1;
      wm->FeederQueue()->wait_and_pop(basket);
      waiting[tid] = 0;
      // Check exit condition: null basket in the queue
      if (!basket) break;
      counter++;
#ifdef __STAT_DEBUG
      sch->GetQueuedStat().RemoveTracks(basket->GetInputTracks());
      sch->GetTransportStat().AddTracks(basket->GetInputTracks());
#endif         
      ntotransport = basket->GetNinput();  // all tracks to be transported 
//      ninput = ntotransport;
      GeantTrack_v &input = basket->GetInputTracks();
      GeantTrack_v &output = basket->GetOutputTracks();
      if (!ntotransport) goto finish;      // input list empty
//      Printf("======= BASKET %p with %d tracks counter=%d =======", basket, ntotransport, counter.load());
//      basket->Print();
//      Printf("==========================================");
//      propagator->fTracksPerBasket[tid] = ntotransport;
      td->fVolume = basket->GetVolume();
      
      // Record tracks
//      ninput = ntotransport;
      if (basket->GetNoutput()) {
         Printf("Ouch: noutput=%d counter=%d", basket->GetNoutput(), counter.load());
      } 
//      if (counter==1) input.PrintTracks();  
/*
      for (Int_t itr=0; itr<ntotransport; itr++) {
         iev[itr] = input.fEventV[itr];
         itrack[itr] = input.fParticleV[itr];
         // crt[itr] = input.fPathV[itr];
         // nxt[itr] = input.fNextpathV[itr];
         if (TMath::IsNaN(input.fXdirV[itr])) {
            Printf("Error: track %d has NaN", itr);
         }   
      }
*/
      // Select the discrete physics process for all particles in the basket
      if (propagator->fUsePhysics) propagator->ProposeStep(ntotransport, input, tid);
      // Apply msc for charged tracks
      propagator->ApplyMsc(ntotransport, input, tid);
      
      ncross = 0;
      generation = 0;
      
      while (ntotransport) {
         // Interrupt condition here. Work stealing could be also implemented here...
         generation++;
//         Printf("====== WorkloadManager:");
//         input.PrintTracks();
         // Propagate all remaining tracks
         if (basket->IsMixed()) 
            ncross += input.PropagateTracksSingle(output,tid);
         else   
            ncross += input.PropagateTracks(output,tid);
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
         nextra_at_rest = 0;
         gPropagator->Process()->Eloss(td->fVolume->GetMaterial(), output.GetNtracks(), output, nextra_at_rest, tid);
//         if (nextra_at_rest) Printf("Extra particles: %d", nextra_at_rest);
      }   
/*
      if (propagator->fUsePhysics) {
         for (Int_t iproc=0; iproc<nprocesses; iproc++) {
            if (propagator->Process(iproc)->IsType(PhysicsProcess::kContinuous)) {
               Int_t nafter = 0;
               gPropagator->Process(iproc)->PostStep(td->fVolume->GetMaterial(), output.GetNtracks(), output, nafter, tid);
               // Now we may also have particles killed by energy threshold
            }   
         }      
      }
*/
      // Now we may also have particles killed by energy threshold
      // Do post-step actions on remaining particles
      // Loop all processes to group particles per process

      if (propagator->fUsePhysics) {
         // Discrete processes only
         Int_t nphys = output.SortByStatus(kPhysics);
         if (nphys) {
            // Do post step actions for particles suffering a given process.
            // Surviving particles are added to the output array
            //propagator->Process()->PostStep(td->fVolume->GetMaterial(), nphys, output, ntotnext, tid);
            
            // first: sample target and type of interaction for each primary tracks 
            propagator->Process()->PostStepTypeOfIntrActSampling( 
                                            td->fVolume->GetMaterial(), 
                                            nphys,
                                            output,
                                            tid);

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
            propagator->fVectorPhysicsProcess->PostStepFinalStateSampling(
                                            td->fVolume->GetMaterial(),
                                            nphys,
                                            output,
                                            ntotnext,
                                            tid);
#endif
            // second: sample final states (based on the inf. regarding sampled
            //         target and type of interaction above), insert them into 
            //         the track vector, update primary tracks; 
            propagator->Process()->PostStepFinalStateSampling(
                                            td->fVolume->GetMaterial(),
                                            nphys,
                                            output,
                                            ntotnext,
                                            tid); 

            if (0 /*ntotnext*/) {
               printf("============= Basket: %s\n", basket->GetName());
               output.PrintTracks();
            }   
         }
      }
      gPropagator->fApplication->StepManager(tid, output.GetNtracks(), output);   
/*
      if (propagator->fUsePhysics) {
         // Discrete processes only
         Int_t nphys = output.SortByStatus(kPhysics);
         if (nphys) {
            for (Int_t iproc=0; iproc<nprocesses; iproc++) {
               if (propagator->Process(iproc)->IsType(PhysicsProcess::kContinuous)) continue;
//               propagator->SelectTracksForProcess(iproc, ntotransport, particles, ntodo, parttodo);
               // Do post step actions for particles suffering a given process.
               // Surviving particles are added to the output array
               propagator->Process(iproc)->PostStep(td->fVolume->GetMaterial(), nphys, output, ntotnext, tid);
            }
         }
      }
*/
      // Check
      if (basket->GetNinput()) {
         Printf("Ouch: ninput=%d counter=%d", basket->GetNinput(), counter.load());
      }   
//      noutput = basket->GetNoutput();
/*
      for(Int_t itr=0; itr<noutput; itr++) {
         if (TMath::IsNaN(output.fXdirV[itr])) {
            Printf("Error: track %d has NaN", itr);
         }   
      }   
      for (Int_t itr=0; itr<ninput; itr++) {         
         Bool_t found = kFALSE;
         for(Int_t i=0; i<noutput; i++) {
            if (output.fEventV[i] == iev[itr]) {
               if (output.fParticleV[i] == itrack[itr]) {
                  found = true;
                  break;
               }   
            }
         }
         if (!found)  {
            Printf("Track %d of event %d not found, counter=%d", itrack[itr],iev[itr], counter.load());
//            output.PrintTracks();
         }   
      }
*/      

finish:
//      basket->Clear();
//      Printf("======= BASKET(tid=%d): in=%d out=%d =======", tid, ninput, basket->GetNoutput());
#ifdef __STAT_DEBUG
      sch->GetTransportStat().RemoveTracks(basket->GetOutputTracks());
#endif         
//      ninjected = 
      sch->AddTracks(basket, ntot, nnew, nkilled);
//      Printf("thread %d: injected %d baskets", tid, ninjected);
//      wm->TransportedQueue()->push(basket);
      wm->SetWorkDone(true);
      cond_sch->notify_one();
      basket->Recycle();
   }
   wm->DoneQueue()->push(0);
   Printf("=== Thread %d: exiting ===", tid);
   return 0;
}

//______________________________________________________________________________
void *WorkloadManager::TransportTracksCoprocessor(void *arg)
{
// Thread propagating all tracks from a basket.
//      char slist[256];
//      TString sslist;
//   const Int_t max_idle = 1;
//   Int_t indmin, indmax;
   static Int_t counter=0;
   //Int_t ntotnext, ncross;
   Int_t ntotransport;
   Int_t generation = 0;
   GeantBasket *basket = 0;
   Int_t tid = TGeoManager::ThreadId();
   GeantPropagator *propagator = GeantPropagator::Instance();
   GeantThreadData *td = propagator->fThreadData[tid];
   WorkloadManager *wm = WorkloadManager::Instance();
   GeantScheduler *sch = wm->GetScheduler();
   Int_t *waiting = wm->GetWaiting();
   //Int_t nprocesses = propagator->fNprocesses;
   // Int_t ninput;
   Int_t noutput;
//   Bool_t useDebug = propagator->fUseDebug;
//   Printf("(%d) WORKER started", tid);
   // Create navigator if none serving this thread.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) nav = gGeoManager->AddNavigator();
   waiting[tid] = 1;
   // Int_t iev[500], itrack[500];
   // TGeoBranchArray *crt[500], *nxt[500];

   TaskBroker *broker = reinterpret_cast<TaskBroker*>(arg);
   while (1) {

      // ::Info("GPU","Waiting (1) for next available stream.");
      TaskBroker::Stream stream = broker->GetNextStream();
      if (!stream) break;
      
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
      //::Info("GPU","Waiting for next available basket.");
      wm->FeederQueue()->wait_and_pop(basket);
      waiting[tid] = 0;
      // Check exit condition: null basket in the queue
      if (!basket) break;
      counter++;

      if (!stream) {
         ::Info("GPU","Waiting (2) for next available stream.");
         stream = broker->GetNextStream();
      }
      // lastToClear = kFALSE;
      if (!basket) {
         if (0 != broker->launchTask(/* wait= */ true)) {
            // We ran something, new basket might be available.
            continue;
         } else {
            break;
         }
      }

      ntotransport = basket->GetNinput();  // all tracks to be transported 
      //ninput = ntotransport;
      GeantTrack_v &input = basket->GetInputTracks();
      GeantTrack_v &output = basket->GetOutputTracks();
      if (!ntotransport) goto finish;      // input list empty
//      Printf("======= BASKET %p with %d tracks counter=%d =======", basket, ntotransport, counter);
//      basket->Print();
//      Printf("==========================================");
//      propagator->fTracksPerBasket[tid] = ntotransport;
      td->fVolume = basket->GetVolume();
      
      // Record tracks
      //ninput = ntotransport;
      if (basket->GetNoutput()) {
         Printf("Ouch (coproc): noutput=%d counter=%d", basket->GetNoutput(), counter);
      } 
//      if (counter==1) input.PrintTracks();  
      for (Int_t itr=0; itr<ntotransport; itr++) {
         // iev[itr] = input.fEventV[itr];
         // itrack[itr] = input.fParticleV[itr];
         // crt[itr] = input.fPathV[itr];
         // nxt[itr] = input.fNextpathV[itr];
         if (TMath::IsNaN(input.fXdirV[itr])) {
            Printf("Error: track %d has NaN", itr);
         }   
      }
      // Select the discrete physics process for all particles in the basket
      //if (propagator->fUsePhysics) propagator->ProposeStep(ntotransport, input, tid);

      //ncross = 0;
      generation = 0;
      
      while (ntotransport) {
         // Interrupt condition here. Work stealing could be also implemented here...
         generation++;
         // Propagate all remaining tracks
         // NOTE: need to deal with propagator->fUsePhysics
         broker->runTask(tid, *basket); // ntotransport, basket_sch->GetNumber(), gPropagator->fTracks, particles);
         ntotransport = 0;
         //ncross += input.PropagateTracks(output);
         //ntotransport = input.GetNtracks();
      }
      // All tracks are now in the output track vector. Possible statuses:
      // kCrossing - particles crossing boundaries
      // kPhysics - particles reaching the point where the discrete physics process 
      //            will happen.
      // kExitingSetup - particles exiting the geometry
      // kKilled - particles that could not advance in geometry after several tries

      // Check
      if (basket->GetNinput()) {
         Printf("Ouch (coproc): ninput=%d noutput=%d counter=%d", basket->GetNinput(), basket->GetNoutput(), counter);
      }   
      noutput = basket->GetNoutput();
      for(Int_t itr=0; itr<noutput; itr++) {
         if (TMath::IsNaN(output.fXdirV[itr])) {
            Printf("Error: track %d has NaN", itr);
         }   
      }   

finish:
//      basket->Clear();
//      Printf("======= BASKET(tid=%d): in=%d out=%d =======", tid, ninput, basket->GetNoutput());
#ifdef __STAT_DEBUG
      sch->GetTransportStat().RemoveTracks(basket->GetOutputTracks());
#endif         
      Int_t ntot = 0;
      Int_t nnew = 0;
      Int_t nkilled = 0;
      /* Int_t ninjected = */ sch->AddTracks(basket, ntot, nnew, nkilled);
      wm->TransportedQueue()->push(basket);
      (void)ntot;
      (void)nnew;
      (void)nkilled;
   }
   wm->DoneQueue()->push(0);
   Printf("=== Coprocessor Thread %d: exiting === Processed %ld", tid, broker->GetTotalWork());
   return 0;
}

//______________________________________________________________________________
void *WorkloadManager::MonitoringThread(void *)
{
// Thread providing basic monitoring for the scheduler.
   const Double_t MByte = 1024.;
   Printf("Started monitoring thread...");
//   GeantPropagator *propagator = GeantPropagator::Instance();
   WorkloadManager *wm = WorkloadManager::Instance();
   Geant::priority_queue<GeantBasket*> *feederQ = wm->FeederQueue();
   Int_t ntotransport;
   ProcInfo_t  procInfo;
   Double_t rss;
   Double_t nmem[100] = {0};
   GeantScheduler *sch = wm->GetScheduler();
   Int_t nvol = sch->GetNvolumes();
   Int_t nthreads = wm->GetNthreads();
   Int_t *nworking = new Int_t[nthreads+1];
   memset(nworking, 0, (nthreads+1)*sizeof(Int_t));
   Int_t *waiting = wm->GetWaiting();
   GeantBasketMgr **bmgr = sch->GetBasketManagers();
   
   TCanvas *cmon = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("cscheduler");
   cmon->Divide(2,2);
   TH1I *hqueue = new TH1I("hqueue","Work queue load", 100,0,100);
   hqueue->SetFillColor(kRed);
   hqueue->SetLineColor(0);
   hqueue->SetStats(kFALSE);
   Int_t nqueue[100] = {0};
   TH1F *hmem = new TH1F("hmem","Resident memory [MB]", 100,0,100);
//   hmem->SetFillColor(kMagenta);
   hmem->SetLineColor(kMagenta);
   hmem->SetStats(kFALSE);
   TH1I *hbaskets = new TH1I("hbaskets","Baskets per volume", nvol,0,nvol);
   hbaskets->SetFillColor(kBlue);
   hbaskets->SetLineColor(0);
   hbaskets->SetStats(kFALSE);
   TH1I *hbused = new TH1I("hbused","Baskets per volume", nvol,0,nvol);
   hbused->SetFillColor(kRed);
   hbused->SetFillStyle(3001);
   hbused->SetLineColor(0);
   hbused->SetStats(kFALSE);
   TH1F *hconcurrency = new TH1F("hconcurrency","Concurrency plot", nthreads+1,0,nthreads+1);
   hconcurrency->GetYaxis()->SetRangeUser(0,1);
   hconcurrency->GetXaxis()->SetNdivisions(nthreads+1,kTRUE);
   hconcurrency->GetYaxis()->SetNdivisions(10,kTRUE);
   hconcurrency->SetFillColor(kGreen);
   hconcurrency->SetLineColor(0);
   hconcurrency->SetStats(kFALSE);
   TH1F *hconcavg = new TH1F("hconcavg","Concurrency plot", nthreads+1,0,nthreads+1);
   hconcavg->SetFillColor(kRed);
   hconcavg->SetFillStyle(3001);
   hconcavg->SetLineColor(0);
   hconcavg->SetStats(kFALSE);
   cmon->cd(1);
   hqueue->Draw();
   cmon->cd(2);
   hmem->Draw();
   cmon->cd(3);
   hbaskets->Draw();
   hbused->Draw("SAME");
   cmon->cd(4);
   hconcurrency->Draw();
   hconcavg->Draw("SAME");
   cmon->Update();
   Double_t stamp = 0.;
   Int_t i, j, bin;
   while (!wm->IsStopped()) { // exit condition here
      i = Int_t(stamp);
      gSystem->Sleep(50); // millisec
      gSystem->GetProcInfo(&procInfo);
      rss = procInfo.fMemResident/MByte;
      ntotransport = feederQ->size_async();
      if (stamp>100) {
         memmove(nqueue, &nqueue[1], 99*sizeof(Int_t));
         memmove(nmem, &nmem[1], 99*sizeof(Double_t));
         nqueue[99] = ntotransport;
         nmem[99] = rss;
         hqueue->GetXaxis()->Set(100,stamp-100, stamp);
         hmem->GetXaxis()->Set(100,stamp-100, stamp);
         for (j=0; j<100; j++) {
            hqueue->SetBinContent(j+1,nqueue[j]);
            hmem->SetBinContent(j+1,nmem[j]);
         }   
      } else {
         nqueue[i] = ntotransport;
         nmem[i] = rss;
         hqueue->SetBinContent(i+1,nqueue[i]);
         hmem->SetBinContent(i+1,nmem[i]);
      } 
      for (j=0; j<nvol; j++) {
         bin = hbaskets->GetXaxis()->FindBin(bmgr[j]->GetVolume()->GetName());
         hbaskets->SetBinContent(bin, bmgr[j]->GetNbaskets());
         hbused->SetBinContent(bin, bmgr[j]->GetNused());
      }
      for (j=0; j<nthreads+1; j++) {
         nworking[j] += 1-waiting[j];
         hconcurrency->SetBinContent(j+1, 1-waiting[j]);
         hconcavg->SetBinContent(j+1, nworking[j]/(stamp+1));
      }   
      cmon->cd(1);
      hqueue->Draw();
      cmon->cd(2);
      hmem->Draw();
      cmon->cd(3);
      hbaskets->LabelsOption(">");
      hbaskets->Draw();
      hbused->Draw("SAME");
      cmon->cd(4);
      hconcurrency->Draw();
      hconcavg->Draw("SAME");
      cmon->Modified();
      cmon->Update();
      stamp += 1;
   }
   delete [] nworking;
   // Sleep a bit to let the graphics finish
   gSystem->Sleep(10); // millisec

   return 0;
}
