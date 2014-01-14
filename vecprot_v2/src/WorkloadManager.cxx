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
#include "TGeoBranchArray.h"
#include "GeantBasket.h"
#include "GeantOutput.h"
#include "GeantThreadData.h"
#include "PhysicsProcess.h"
#include "GeantScheduler.h"
#include "GeantEvent.h"

ClassImp(WorkloadManager)

WorkloadManager *WorkloadManager::fgInstance = 0;

//______________________________________________________________________________
WorkloadManager::WorkloadManager(Int_t nthreads)
{
// Private constructor.
   fNthreads = nthreads;
   fNbaskets = 0;
   fBasketGeneration = 0;
   fNbasketgen = 0;
   fNidle = nthreads;
   fNminThreshold = 10;
   fNqueued = 0;
   fBtogo = 0;
   fStarted = kFALSE;
   fFeederQ = new dcqueue<GeantBasket>();
   fTransportedQ = new dcqueue<GeantBasket>();
   fDoneQ = new dcqueue<GeantBasket>();
   fgInstance = this;
   fListThreads = 0;
   fFlushed = kFALSE;
   fFilling = kFALSE;
   fScheduler = new GeantScheduler();
}

//______________________________________________________________________________
WorkloadManager::~WorkloadManager()
{
// Destructor.
   delete fFeederQ;
   delete fTransportedQ;
   delete fDoneQ;
   fgInstance = 0;
}

//______________________________________________________________________________
void WorkloadManager::CreateBaskets()
{
// Create the array of baskets
   fScheduler->CreateBaskets();
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
void WorkloadManager::Print(Option_t *option) const
{
//
}   

//______________________________________________________________________________
void WorkloadManager::StartThreads()
{
// Start the threads
   fStarted = kTRUE;
   if (fListThreads) return;
   fListThreads = new TList();
   fListThreads->SetOwner();
   for (Int_t ith=0; ith<fNthreads; ith++) {
      TThread *t = new TThread(WorkloadManager::TransportTracks);
      fListThreads->Add(t);
      t->Run();
   }
//   gSystem->Sleep(1000);
   TThread *t = new TThread(WorkloadManager::MainScheduler);
   fListThreads->Add(t);
   t->Run();
}   

//______________________________________________________________________________
void WorkloadManager::JoinThreads()
{
// 
   for (Int_t ith=0; ith<fNthreads; ith++) fFeederQ->push(0);
   for (Int_t ith=0; ith<fNthreads; ith++) ((TThread*)fListThreads->At(ith))->Join();
   // Join garbage collector
   ((TThread*)fListThreads->At(fNthreads))->Join();
}
   
//______________________________________________________________________________
void WorkloadManager::WaitWorkers()
{
// Waiting point for the main thread until work gets done.
   fFilling = kFALSE;
   Int_t ntowait = fNthreads;
   while(ntowait) {
      fDoneQ->wait_and_pop();
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
   
   Bool_t graphics = propagator->fUseGraphics;
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
   dcqueue<GeantBasket> *feederQ = wm->FeederQueue();
   dcqueue<GeantBasket> *transportedQ = wm->TransportedQueue();
   GeantScheduler *sch = wm->GetScheduler();
   // Number of baskets in the queue to transport
   Int_t ntotransport = 0;
//   Int_t ntracksperbasket = propagator->fNperBasket;
   // Feeder threshold
   Int_t min_feeder = TMath::Max(20,3*nworkers); // setter/getter here ?
   Int_t max_feeder = 300; // ???
   // Number of tracks in the current basket
//   Int_t ntracksb = 0;
   // Number of priority baskets
   Int_t npriority;
   TH1F *hnb=0, *hworkers=0, *htracks=0;
   TCanvas *c1 = 0;
   Int_t lastphase = -1;
   Int_t crtphase  = 0;
   
   if (graphics) {
      hnb = new TH1F("hnb","number of baskets in the transport queue; iteration#",200000,0,200000);
      hnb->SetFillColor(kRed);
      hnb->SetLineColor(0);
      hworkers = new TH1F("hworkers","number of active workers; iteration#",200000,0,200000);
      hworkers->SetLineColor(0);
      hworkers->SetFillColor(kBlue);
      htracks = new TH1F("htracks","number of tracks/basket; iteration#",200000,0,200000);
      htracks->SetLineColor(0);
      htracks->SetFillColor(kGreen);
   }      
   Int_t niter = -1;
   UInt_t npop = 0;
   Double_t nperbasket;
   Int_t ninjected = 0;
   Int_t nwaiting = 0;
   Int_t nnew = 0, ntot=0, nkilled=0;
   GeantBasket *output;
   GeantBasket **carray = new GeantBasket*[500];
   while ((output = (GeantBasket*)transportedQ->wait_and_pop_max(500,npop,carray))) {
      // Monitor the queues while there are tracks to transport
      niter++;
      ninjected = 0;
      nnew = 0;
      ntot = 0;
      nkilled = 0;
      // Check first the queue of collectors
//      ncollectors = npop+collector_queue->size_async();
//      Printf("Popped %d collections, %d still in the queue", npop, ncollectors-npop);
//      nempty = empty_queue->size_async();
      // Process popped collections and flush their tracks
      for (UInt_t iout=0; iout<npop; iout++) {
         output = carray[iout];
//         collector->Print();
//         Printf("= collector has %d tracks", collector->GetNtracks());
         ninjected += sch->AddTracks(output, ntot, nnew, nkilled);
//         Printf("=== injected %d baskets", ninjected);
         // Recycle basket
	      output->Recycle();	 
      }
      // If there were events to be dumped, check their status here
      ntotransport = feederQ->size_async();
      Printf("#%d: feeder=%p Processed %d baskets (%d tracks, %d new, %d killed)-> injected %d. QS=%d", niter, feederQ, npop, ntot, nnew, nkilled, ninjected, ntotransport);
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
//            for (Int_t itrack=0; itrack<ntracks; itrack++) {
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
                  Printf("Flushed %d priority baskets, resetting countdown", npriority);
               } else {
                  Printf("Countdown is %d", feederQ->get_countdown());
               }
            } else {
               npriority = sch->FlushPriorityBaskets();
               ninjected = npriority;
               Printf("Flushed %d priority baskets", npriority);
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
      nwaiting = propagator->GetNwaiting();
//      Printf("picked=%d ncoll=%d  ninjected=%d ntotransport=%d",npop, ncollectors,ninjected,ntotransport);
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
           continue;
        }
        nwaiting = propagator->GetNwaiting();
//        if (!ninjected && !ntotransport && !nnew && nwaiting==nworkers) break;
//        direct_feed = (ntotransport==0)?kTRUE:kFALSE;
      }
      ntotransport = feederQ->size_async();
      if (ntotransport==0) {
         Printf("Garbage collection");
         sch->GarbageCollect();
         if (countdown) feederQ->set_countdown(0);
      }
      nperbasket = 0;
//      for (Int_t tid=0; tid<propagator->fNthreads; tid++) nperbasket += propagator->fTracksPerBasket[tid];
//      nperbasket /= propagator->fNthreads;
      if (graphics) {
         if ((niter%10==0) && (niter/10<200000)) {
            hnb->Fill(niter/10, ntotransport);
            hworkers->Fill(niter/10, nworkers-nwaiting);
            htracks->Fill(niter/10, nperbasket);
         }   
      }   
   }
      
   Printf("=== Scheduler: stopping threads === niter =%d\n", niter);
   if (graphics) {
      c1 = new TCanvas("c2","c2",1200,1200);
      c1->Divide(1,3);
      c1->cd(1);
      hnb->SetStats(kFALSE);
      hnb->GetXaxis()->SetRangeUser(0,niter/10);
      hnb->Draw();
      c1->cd(2);
      hworkers->SetStats(kFALSE);
      hworkers->GetXaxis()->SetRangeUser(0,niter/10);
      hworkers->Draw();
      c1->cd(3);
      htracks->SetStats(kFALSE);
      htracks->GetXaxis()->SetRangeUser(0,niter/10);
      htracks->Draw();
//      pad1->Modified();
//      pad2->Modified();
//      pad3->Modified();
//      c1->Update();
   }   
   for (Int_t i=0; i<propagator->fNthreads; i++) wm->FeederQueue()->push(0);
   wm->TransportedQueue()->push(0);
   Printf("=== Scheduler: exiting ===");
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
   static Int_t counter=0;
   Int_t ntotnext, ncross, ntotransport;
   Int_t generation = 0;
   GeantBasket *basket = 0;
   Int_t tid = TGeoManager::ThreadId();
   GeantPropagator *propagator = GeantPropagator::Instance();
   GeantThreadData *td = propagator->fThreadData[tid];
   WorkloadManager *wm = WorkloadManager::Instance();
   GeantScheduler *sch = wm->GetScheduler();
   Int_t nprocesses = propagator->fNprocesses;
   Int_t ninput, noutput;
//   Bool_t useDebug = propagator->fUseDebug;
//   Printf("(%d) WORKER started", tid);
   // Create navigator if none serving this thread.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) nav = gGeoManager->AddNavigator();
   propagator->fWaiting[tid] = 1;
   Int_t iev[500], itrack[500];
   TGeoBranchArray *crt[500], *nxt[500];
   while (1) {
      propagator->fWaiting[tid] = 1;
      basket = wm->FeederQueue()->wait_and_pop();
      propagator->fWaiting[tid] = 0;
      // Check exit condition: null basket in the queue
      if (!basket) break;
      counter++;
#ifdef __STAT_DEBUG
      sch->GetQueuedStat().RemoveTracks(basket->GetInputTracks());
      sch->GetTransportStat().AddTracks(basket->GetInputTracks());
#endif         
      ntotransport = basket->GetNinput();  // all tracks to be transported 
      ninput = ntotransport;
      GeantTrack_v &input = basket->GetInputTracks();
      GeantTrack_v &output = basket->GetOutputTracks();
      if (!ntotransport) goto finish;      // input list empty
//      Printf("======= BASKET %p with %d tracks counter=%d =======", basket, ntotransport, counter);
//      basket->Print();
//      Printf("==========================================");
//      propagator->fTracksPerBasket[tid] = ntotransport;
      td->fVolume = basket->GetVolume();
      
      // Record tracks
      ninput = ntotransport;
      if (basket->GetNoutput()) {
         Printf("Ouch: noutput=%d counter=%d", basket->GetNoutput(), counter);
      } 
//      if (counter==1) input.PrintTracks();  
      for (Int_t itr=0; itr<ntotransport; itr++) {
         iev[itr] = input.fEventV[itr];
         itrack[itr] = input.fParticleV[itr];
         crt[itr] = input.fPathV[itr];
         nxt[itr] = input.fNextpathV[itr];
         if (TMath::IsNaN(input.fXdirV[itr])) {
            Printf("Error: track %d has NaN", itr);
         }   
      }
      // Select the discrete physics process for all particles in the basket
      if (propagator->fUsePhysics) propagator->PhysicsSelect(ntotransport, input, tid);
      
      ncross = 0;
      generation = 0;
      
      while (ntotransport) {
         // Interrupt condition here. Work stealing could be also implemented here...
         generation++;
         // Propagate all remaining tracks
         ncross += input.PropagateTracks(output);
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
         for (Int_t iproc=0; iproc<nprocesses; iproc++) {
            if (propagator->Process(iproc)->IsType(PhysicsProcess::kContinuous)) {
               Int_t nafter = 0;
               gPropagator->Process(iproc)->PostStep(td->fVolume, output.GetNtracks(), output, nafter, tid);
               // Now we may also have particles killed by energy threshold
            }   
         }      
      }
      // Now we may also have particles killed by energy threshold
      // Do post-step actions on remaining particles
      // Loop all processes to group particles per process
      if (propagator->fUsePhysics) {
         // Discrete processes only
         Int_t nphys = output.SortByStatus(kPhysics);
         if (nphys) {
            for (Int_t iproc=0; iproc<nprocesses; iproc++) {
               if (propagator->Process(iproc)->IsType(PhysicsProcess::kContinuous)) continue;
//               propagator->SelectTracksForProcess(iproc, ntotransport, particles, ntodo, parttodo);
               // Do post step actions for particles suffering a given process.
               // Surviving particles are added to the output array
               propagator->Process(iproc)->PostStep(td->fVolume, nphys, output, ntotnext, tid);
            }
         }
      }
      // Check
      if (basket->GetNinput()) {
         Printf("Ouch: ninput=%d counter=%d", basket->GetNoutput(), counter);
      }   
      noutput = basket->GetNoutput();
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
            Printf("Track %d of event %d not found, counter=%d", itrack[itr],iev[itr], counter);
//            output.PrintTracks();
         }   
      }

finish:
//      basket->Clear();
//      Printf("======= BASKET(tid=%d): in=%d out=%d =======", tid, ninput, basket->GetNoutput());
#ifdef __STAT_DEBUG
      sch->GetTransportStat().RemoveTracks(basket->GetOutputTracks());
#endif         
      wm->TransportedQueue()->push(basket);
   }
   wm->DoneQueue()->push(0);
   Printf("=== Thread %d: exiting ===", tid);
   return 0;
}
