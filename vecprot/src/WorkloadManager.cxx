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
#include "GeantVolumeBasket.h"
#include "GeantOutput.h"
#include "GeantThreadData.h"
#include "PhysicsProcess.h"

#include "TaskBroker.h"

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
   answer_queue = new concurrent_queue();
   feeder_queue = new concurrent_queue(true);
   empty_queue = new concurrent_queue;
   collector_queue = new concurrent_queue();
   collector_empty_queue = new concurrent_queue();
   for (Int_t i=0; i<nthreads; i++) collector_empty_queue->push(new GeantTrackCollection(100));
   fStarted = kFALSE;
   fgInstance = this;
   fCurrentBasket = new GeantVolumeBasket*[nthreads];
   fListThreads = 0;
   fBasketArray = 0;
   fFlushed = kFALSE;
   fFilling = kFALSE;
   fBroker = 0;
}

//______________________________________________________________________________
WorkloadManager::~WorkloadManager()
{
   // Destructor.
   delete answer_queue;
   delete feeder_queue;
   delete empty_queue;
   delete [] fCurrentBasket;
   if (fBasketArray) {
      for (Int_t i=0; i<fNbaskets; i++) delete fBasketArray[i];
      delete [] fBasketArray;
   }
   fgInstance = 0;
}

//______________________________________________________________________________
void WorkloadManager::CreateBaskets(Int_t nvolumes)
{
   // Create the array of baskets
   if (fBasketArray) return;
   fBasketArray = new GeantVolumeBasket*[nvolumes];
   TIter next(gGeoManager->GetListOfVolumes());
   TGeoVolume *vol;
   GeantVolumeBasket *basket;
   Int_t icrt = 0;
   while ((vol=(TGeoVolume*)next())) {
      basket = new GeantVolumeBasket(vol, icrt++);
      vol->SetField(basket);
      AddBasket(basket);
   }
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
   feeder_queue->Print();
}

//______________________________________________________________________________
void WorkloadManager::AddEmptyBaskets(Int_t nb)
{
   // Add empty baskets in the queue.
   for (Int_t i=0; i<nb; i++) empty_queue->push(new GeantBasket(10));
   Printf("Added %d empty baskets to the queue", nb);
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
      Printf("Running with a coprocessor broker.");
      TThread *t = new TThread(WorkloadManager::TransportTracksCoprocessor,fBroker);
      fListThreads->Add(t);
      t->Run();
      ++ith;
   }
   for (; ith<fNthreads; ith++) {
      TThread *t = new TThread(WorkloadManager::TransportTracks);
      fListThreads->Add(t);
      t->Run();
   }
   gSystem->Sleep(1000);
   TThread *t = new TThread(WorkloadManager::MainScheduler);
   fListThreads->Add(t);
   t->Run();
}

//______________________________________________________________________________
void WorkloadManager::JoinThreads()
{
   //
   for (Int_t ith=0; ith<fNthreads; ith++) feeder_queue->push(0);
   for (Int_t ith=0; ith<fNthreads; ith++) ((TThread*)fListThreads->At(ith))->Join();
   // Join garbage collector
   ((TThread*)fListThreads->At(fNthreads))->Join();
}

//______________________________________________________________________________
void WorkloadManager::WaitWorkers()
{
   // Waiting point for the main thread until work gets done.
   fFilling = kFALSE;
   Int_t ntowait = fNthreads+1;
   while(ntowait) {
      answer_queue->wait_and_pop();
      ntowait--;
      //      Printf("Worker %d finished", ipop);
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
   WorkloadManager *wm = WorkloadManager::Instance();
   concurrent_queue *feeder_queue = wm->FeederQueue();
   concurrent_queue *empty_queue = wm->EmptyQueue();
   concurrent_queue *collector_queue = wm->CollectorQueue();
   concurrent_queue *collector_empty_queue = wm->CollectorEmptyQueue();
   // Number of baskets in the queue to transport
   Int_t ntotransport = 0;
   // Number of collectors available
   Int_t ncollectors = 0;
   // Number of empty baskets available
   Int_t nempty = 0;
   //   Int_t ntracksperbasket = propagator->fNperBasket;
   // Feeder threshold
   Int_t min_feeder = TMath::Max(50,2*propagator->fNthreads); // setter/getter here ?
   // Number of tracks in the current basket
   //   Int_t ntracksb = 0;
   Int_t nbaskets = wm->GetNbaskets();
   // Number of priority baskets
   Int_t npriority;
   // Main scheduler
   GeantMainScheduler *main_sch = new GeantMainScheduler(nbaskets);
   TH1F *hnb=0, *hworkers=0, *htracks=0;
   TCanvas *c1 = 0;
   TPad *pad1=0, *pad2=0, *pad3=0;
   Int_t lastphase = -1;
   Int_t crtphase  = 0;
   
   if (graphics) {
      hnb = new TH1F("hnb","number of baskets in the transport queue; iteration#",200000,0,200000);
      hnb->SetFillColor(kRed);
      hworkers = new TH1F("hworkers","number of active workers; iteration#",200000,0,200000);
      hworkers->SetFillColor(kBlue);
      htracks = new TH1F("htracks","number of tracks/basket; iteration#",200000,0,200000);
      htracks->SetFillColor(kGreen);
   }
   Int_t niter = -1;
   UInt_t npop = 0;
   Double_t nperbasket;
   Int_t ninjected = 0;
   Int_t nwaiting = 0;
   Int_t nnew = 0;
   GeantTrackCollection *collector;
   TObject **carray = new TObject*[500];
   while ((collector = (GeantTrackCollection*)collector_queue->wait_and_pop_max(500,npop,carray))) {
      // Monitor the queues while there are tracks to transport
      niter++;
      ninjected = 0;
      nnew = 0;
      // Check first the queue of collectors
      ncollectors = npop+collector_queue->size_async();
      //      Printf("Popped %d collections, %d still in the queue", npop, ncollectors-npop);
      nempty = empty_queue->size_async();
      // Process popped collections and flush their tracks
      for (UInt_t icoll=0; icoll<npop; icoll++) {
         collector = (GeantTrackCollection*)carray[icoll];
         //         collector->Print();
         //         Printf("= collector has %d tracks", collector->GetNtracks());
         ninjected += collector->FlushTracks(main_sch);
         //         Printf("=== injected %d baskets", ninjected);
         collector_empty_queue->push(collector);
      }
      // If there were events to be dumped, check their status here
      ntotransport = feeder_queue->size_async();
      // Print info about current phase
      if (crtphase!=lastphase) {
         lastphase=crtphase;
         switch (crtphase) {
            case 0:
               Printf("============         Phase 0:      ============");
               Printf("   Propagating initial population of %d baskets.", ntotransport);
               break;
            case 1:
               Printf("============         Phase 1:      ============");
               Printf("   Dumping ranges of events");
               break;
            case 2:
               Printf("============         Phase 2:      ============");
               Printf("   Importing new events");
               break;
            default:
               break;
         }
      }
      
      // Check and mark finished events
      for (Int_t ievt=0; ievt<nbuffered; ievt++) {
         GeantEvent *evt = propagator->fEvents[ievt];
         if (finished.TestBitNumber(evt->event)) continue;
         if (evt->Transported()) {
            // Digitizer (delete for now)
            Int_t ntracks = propagator->fNtracks[ievt];
            Printf("= digitizing event %d with %d tracks", evt->event, ntracks);
            //            for (Int_t itrack=0; itrack<ntracks; itrack++) {
            //               delete propagator->fTracks[maxperevent*ievt+itrack];
            //               propagator->fTracks[maxperevent*ievt+itrack] = 0;
            //            }
            finished.SetBitNumber(evt->event);
            if (last_event<max_events) {
               Printf("=> Importing event %d", last_event);
               propagator->ImportTracks(1,propagator->fNaverage,last_event,ievt);
               last_event++;
               nnew++;
            }
         }
      }
      if (finished.FirstNullBit() >= (UInt_t)max_events) {
         // Printf("All Events are finished.");
         break;
      }
      // else {
      //    finished.Print("all");
      // }
      
      // In case events were transported with priority, check if they finished
      if (prioritize) {
         first_not_transported = finished.FirstNullBit();
         if (first_not_transported > dumped_event+4) {
            // Priority events digitized, exit prioritized regime
            Printf("= stopped prioritizing");
            prioritize = kFALSE;
            main_sch->SetPriorityRange(-1, -1);
         } else {
            // Flush priority baskets
            npriority = main_sch->GetNpriority();
            if (npriority) ninjected += main_sch->FlushPriorityBaskets();
         }
      }
      
      ntotransport = feeder_queue->size_async();
      nwaiting = propagator->GetNwaiting();
      //      Printf("picked=%d ncoll=%d  ninjected=%d ntotransport=%d",npop, ncollectors,ninjected,ntotransport);
      if (ntotransport<min_feeder) {
         // Transport queue below the threshold
         if (crtphase<1) crtphase = 1;
         // In case no new events were injected and we are not in a the priority regime
         // and below lowest watermark in this iteration, make a garbage collection
         //        if (!nnew && !prioritize && ntotransport<nworkers) {
         ninjected += main_sch->FlushBaskets();
         //           Printf("Garbage collection injected %d baskets", ninjected);
         //        }
         if (!prioritize && last_event<max_events && nnew==0) {
            // Start prioritized regime
            dumped_event = finished.FirstNullBit();
            Printf("Prioritizing events %d to %d", dumped_event,dumped_event+4);
            main_sch->SetPriorityRange(dumped_event, dumped_event+4);
            prioritize = kTRUE;
            continue;
         }
         nwaiting = propagator->GetNwaiting();
         //        if (!ninjected && !ntotransport && !nnew && nwaiting==nworkers) break;
         //        direct_feed = (ntotransport==0)?kTRUE:kFALSE;
      }
      nperbasket = 0;
      for (Int_t tid=0; tid<propagator->fNthreads; tid++) nperbasket += propagator->fTracksPerBasket[tid];
      nperbasket /= propagator->fNthreads;
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
      pad1 = (TPad*)c1->cd(1);
      hnb->SetStats(kFALSE);
      hnb->GetXaxis()->SetRangeUser(0,niter/10);
      hnb->Draw();
      pad2 = (TPad*)c1->cd(2);
      hworkers->SetStats(kFALSE);
      hworkers->GetXaxis()->SetRangeUser(0,niter/10);
      hworkers->Draw();
      pad3 = (TPad*)c1->cd(3);
      htracks->SetStats(kFALSE);
      htracks->GetXaxis()->SetRangeUser(0,niter/10);
      htracks->Draw();
      //      pad1->Modified();
      //      pad2->Modified();
      //      pad3->Modified();
      //      c1->Update();
   }
   for (Int_t i=0; i<propagator->fNthreads+1; i++) wm->FeederQueue()->push(0);
   wm->AnswerQueue()->push(0);
   Printf("=== Scheduler: exiting ===");
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
   Int_t ntotnext, ntmp, ntodo, ncross, cputime, ntotransport;
   GeantTrack *track = 0;
   Int_t itrack;
   Int_t generation = 0;
   Bool_t lastToClear = kFALSE;
   GeantVolumeBasket *basket_sch = 0;
   GeantBasket *basket = 0;
   Int_t tid = TGeoManager::ThreadId();
   GeantPropagator *propagator = GeantPropagator::Instance();
   GeantThreadData *td = propagator->fThreadData[tid];
   WorkloadManager *wm = WorkloadManager::Instance();
   Int_t *particles = td->fPartInd->GetArray();
   Int_t *partnext  = td->fPartNext->GetArray();
   Int_t *parttodo  = td->fPartTodo->GetArray();
   Int_t *partcross = td->fPartCross->GetArray();
   Int_t nprocesses = propagator->fNprocesses;
   //   Bool_t useDebug = propagator->fUseDebug;
   GeantTrack **tracks = propagator->fTracks;
   TGeoBranchArray *path = 0;
   //   Printf("(%d) WORKER started", tid);
   // Create navigator if none serving this thread.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) nav = gGeoManager->AddNavigator();
   propagator->fWaiting[tid] = 1;
   TaskBroker *broker = reinterpret_cast<TaskBroker*>(arg);
   while (1) {
      // fprintf(stderr,"DEBUG2: check slot 55: event# %d %p\n",tracks[55]->event,tracks[55]);

      propagator->fWaiting[tid] = 1;
      ::Info("GPU","Grabbing stuff");
      basket = (GeantBasket*)wm->FeederQueue()->wait_and_pop();
      // gSystem->Sleep(3000);
      propagator->fWaiting[tid] = 0;
      lastToClear = kFALSE;
      if (!basket) break;
      ntotransport = basket->GetNtracks();  // all tracks to be transported
      if (!ntotransport) goto finish;
      //      Printf("======= BASKET %p taken by thread #%d =======", basket, tid);
      //      basket->Print();
      //      Printf("==========================================");
      propagator->fTracksPerBasket[tid] = ntotransport;
      path = tracks[basket->GetTracks()[0]]->path;
      td->fVolume = path->GetCurrentNode()->GetVolume();
      basket_sch = (GeantVolumeBasket*)td->fVolume->GetField();
      wm->SetCurrentBasket(tid,basket_sch);
      Printf("(%d - GPU) ================= BASKET of %s (%d): %d tracks", tid, basket_sch->GetName(), basket_sch->GetNumber(), ntotransport);
      memcpy(particles, basket->GetTracks(), ntotransport*sizeof(Int_t));
      // for(int ti = 0; ti < basket->GetNtracks(); ++ti) {
      //    fprintf(stderr,"DEBUG7: %d points to %d\n", ti, basket->GetTracks()[ti]);
      // }
      //      PrintParticles(particles, ntotransport, tid);
      //      sprintf(slist,"Thread #%d transporting %d tracks in basket %s: ", tid, ntotransport, basket->GetName());
      //      sslist = slist;
      //      for (Int_t ip=0; ip<ntotransport; ip++) {sprintf(slist,"%d ",particles[ip]), sslist += slist;}
      //      Printf("%s", sslist.Data());
      // for (Int_t ip=0; ip<ntotransport; ip++) {
      //    int ti = particles[ip];
      //    fprintf(stderr,"DEBUG track %d:%d is %d \n",tracks[ti]->event,ti,tracks[ti]->IsAlive()); 
      // }
      ntotnext = 0;
      ntmp = 0;
      ntodo = 0;
      ncross = 0;
      cputime = 0.;
      generation = 0;
      track = 0;
      while (ntotransport) {
         generation++;
         // Loop all tracks to generate physics/geometry steps
         // Physics step
         
         broker->runTask(tid, ntotransport, basket_sch->GetNumber(), gPropagator->fTracks, particles);
         
         // Copy only tracks that survived boundaries (well we will have to think of
         // those too, like passing them to the next volume...)
         memcpy(particles, partnext, ntotnext*sizeof(Int_t));
         
         for (Int_t t = 0; t < ntotransport; ++t) {
            if (particles[t]>=0) {
//               if (gPropagator->fTracks[particles[t]]->path->GetLevel() > 1) {
//                  fprintf(stderr,"DEBUG: for %d level is %d\n",t,gPropagator->fTracks[particles[t]]->path->GetLevel());
//               }
               propagator->StopTrack(gPropagator->fTracks[particles[t]]);
            }
         }
         ntotransport = 0;
         
         // I/O: Dump current generation
         //         Printf("   ### Generation %d:  %d tracks  cputime=%f", generation, ntotransport,cputime);
         if (propagator->fFillTree) {
            cputime = gPropagator->fTimer->CpuTime();
            gPropagator->fOutput->SetStamp(gPropagator->fThreadData[tid]->fVolume->GetNumber(), propagator->fWMgr->GetBasketGeneration(), generation, ntotransport, cputime);
            for (itrack=0; itrack<ntotransport;itrack++) {
               if (particles[itrack]>=0) {
                  track = tracks[particles[itrack]];
                  gPropagator->fOutput->SetTrack(itrack, track);
               }
            }
            gPropagator->fOutTree->Fill();
         }
      }
      
   finish:
      ::Info("GPU","Clearing");
      basket->Clear();
      wm->EmptyQueue()->push(basket);
      propagator->InjectCollection(tid);
   }
   wm->AnswerQueue()->push(0);
   Printf("=== Thread %d: exiting ===", tid);
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
   Int_t ntotnext, ntmp, ntodo, ncross, cputime, ntotransport;
   GeantTrack *track = 0;
   Int_t itrack;
   Int_t generation = 0;
   Bool_t lastToClear = kFALSE;
   GeantVolumeBasket *basket_sch = 0;
   GeantBasket *basket = 0;
   Int_t tid = TGeoManager::ThreadId();
   GeantPropagator *propagator = GeantPropagator::Instance();
   GeantThreadData *td = propagator->fThreadData[tid];
   WorkloadManager *wm = WorkloadManager::Instance();
   Int_t *particles = td->fPartInd->GetArray();
   Int_t *partnext  = td->fPartNext->GetArray();
   Int_t *parttodo  = td->fPartTodo->GetArray();
   Int_t *partcross = td->fPartCross->GetArray();
   Int_t nprocesses = propagator->fNprocesses;
   //   Bool_t useDebug = propagator->fUseDebug;
   GeantTrack **tracks = propagator->fTracks;
   TGeoBranchArray *path = 0;
   //   Printf("(%d) WORKER started", tid);
   // Create navigator if none serving this thread.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) nav = gGeoManager->AddNavigator();
   propagator->fWaiting[tid] = 1;
   while (1) {
      // fprintf(stderr,"DEBUG: check slot  5: event# %d %p\n",tracks[ 5]->event,tracks[ 5]);
      // fprintf(stderr,"DEBUG: check slot 55: event# %d %p\n",tracks[55]->event,tracks[55]);

      propagator->fWaiting[tid] = 1;
      basket = (GeantBasket*)wm->FeederQueue()->wait_and_pop();
      propagator->fWaiting[tid] = 0;
      lastToClear = kFALSE;
      if (!basket) {
         // Printf("Finishing for %d with %p\n",tid,basket);
         break;
      }
      ntotransport = basket->GetNtracks();  // all tracks to be transported
      if (!ntotransport) {
         // Printf("Finishing for %d with %p\n",tid,basket);
         goto finish;
      }
      //      Printf("======= BASKET %p taken by thread #%d =======", basket, tid);
      //      basket->Print();
      //      Printf("==========================================");
      propagator->fTracksPerBasket[tid] = ntotransport;
      path = tracks[basket->GetTracks()[0]]->path;
      td->fVolume = path->GetCurrentNode()->GetVolume();
      basket_sch = (GeantVolumeBasket*)td->fVolume->GetField();
      wm->SetCurrentBasket(tid,basket_sch);
      Printf("(%d - CPU) ================= BASKET of %s (%d): %d tracks", tid, basket_sch->GetName(), basket_sch->GetNumber(), ntotransport);
      //      Printf("(%d) ================= BASKET of %s: %d tracks", tid, gPropagator->fVolume[tid]->GetName(), ntotransport);

      // for(int ti = 0; ti < basket->GetNtracks(); ++ti) {
      //    fprintf(stderr,"DEBUG6: %d points to %d\n", ti, basket->GetTracks()[ti]);
      // }
      memcpy(particles, basket->GetTracks(), ntotransport*sizeof(Int_t));
      //      PrintParticles(particles, ntotransport, tid);
      //      sprintf(slist,"Thread #%d transporting %d tracks in basket %s: ", tid, ntotransport, basket->GetName());
      //      sslist = slist;
      //      for (Int_t ip=0; ip<ntotransport; ip++) {sprintf(slist,"%d ",particles[ip]), sslist += slist;}
      //      Printf("%s", sslist.Data());
      // for (Int_t ip=0; ip<ntotransport; ip++) {
      //    int ti = particles[ip];
      //    fprintf(stderr,"DEBUG track %d:%d is %d \n",tracks[ti]->event,ti,tracks[ti]->IsAlive()); 
      // }
      ntotnext = 0;
      ntmp = 0;
      ntodo = 0;
      ncross = 0;
      cputime = 0.;
      generation = 0;
      track = 0;
      while (ntotransport) {
         // Interrupt condition here. Work stealing could be also implemented here...
         /*
          if (!wm->IsFilling() && wm->FeederQueue()->empty()) {
          Int_t nworkers = wm->FeederQueue()->assigned_workers();
          if (wm->GetNthreads()-nworkers > max_idle) {
          Printf("Interrupting transport of %d tracks in %s", ntotransport, basket->GetName());
          wm->InterruptBasket(basket, particles, ntotransport, tid);
          }
          }
          */
         generation++;
         // Loop all tracks to generate physics/geometry steps
         // Physics step
         if (propagator->fUsePhysics) propagator->PhysicsSelect(ntotransport, particles, tid);
         // Geometry snext and safety
         basket_sch->ComputeTransportLength(ntotransport, particles);
         // Propagate tracks with physics step and check if they survive boundaries
         // or physics
         ntmp = ntotransport;
         ntodo = ntotransport;
         ntotnext = 0;
         Int_t *ptrParticles = particles;
         // Propagate all tracks alive with physics step.
         // If a boundary is encountered before the physics step, stop the track
         //         if (useDebug)
         //            Printf("(%d) --- propagating %d tracks for volume basket %s", tid, ntodo, basket_sch->GetName());
         while (ntodo) {
            ntodo = 0;
            ncross = 0;
            // Propagate ALL ntmp tracks
            basket_sch->PropagateTracks(ntmp, ptrParticles, ntotnext, partnext, ntodo, parttodo, ncross, partcross);
            //            printf("(%d) ->crossing particles (%d): ",tid, ncross);
            //            for (Int_t ii=0; ii<ncross; ii++) printf("%d ", partcross[ii]);
            //            printf("\n(%d) ->remaining particles (%d): ", tid, ntotnext);
            //            for (Int_t ii=0; ii<ntotnext; ii++) printf("%d ", partnext[ii]);
            //            printf("\n(%d) ->todo particles (%d): ", tid,ntodo);
            //            for (Int_t ii=0; ii<ntodo; ii++) printf("%d ", parttodo[ii]);
            //            printf("\n");
            // Post-step actions by continuous processes for particles reaching boundaries
            if (propagator->fUsePhysics && ncross) {
               for (Int_t iproc=0; iproc<nprocesses; iproc++) {
                  if (propagator->Process(iproc)->IsType(PhysicsProcess::kDiscrete)) continue;
                  Int_t nafter = 0;
                  gPropagator->Process(iproc)->PostStep(td->fVolume, ncross, partcross, nafter, NULL,tid);
                  basket_sch->ResetStep(ncross, partcross);
               }
            }
            ntmp = ntodo;
            ptrParticles = parttodo;
         }
         
         // Copy only tracks that survived boundaries (well we will have to think of
         // those too, like passing them to the next volume...)
         memcpy(particles, partnext, ntotnext*sizeof(Int_t));
         ntotransport = ntotnext;
         
         // Do post-step actions on remaining particles
         ntotnext = 0;
         // Loop all processes to group particles per process
         if (propagator->fUsePhysics && ntotransport) {
            // Apply continuous processes to all particles
            for (Int_t iproc=0; iproc<nprocesses; iproc++) {
               if (propagator->Process(iproc)->IsType(PhysicsProcess::kDiscrete)) continue;
               ntodo = 0;
               gPropagator->Process(iproc)->PostStep(td->fVolume, ntotransport, particles, ntodo, parttodo, tid);
               // Do we have stopped particles ?
               if (ntodo<ntotransport) {
                  memcpy(particles, parttodo, ntodo*sizeof(Int_t));
                  ntotransport = ntodo;
               }
            }
            // Copy al tracks for which step was limited by a continuous process
            // to the next array
            for (Int_t itr=0; itr<ntotransport; itr++) {
               if (tracks[particles[itr]]->process < 0) continue;
               if (propagator->Process(tracks[particles[itr]]->process)->IsType(PhysicsProcess::kContinuous))
                  partnext[ntotnext++] = particles[itr];
            }
            // Discrete processes only
            for (Int_t iproc=0; iproc<nprocesses; iproc++) {
               // Make arrays of particles per process -> ntodo, parttodo
               if (propagator->Process(iproc)->IsType(PhysicsProcess::kContinuous)) continue;
               ntodo = 0;
               propagator->SelectTracksForProcess(iproc, ntotransport, particles, ntodo, parttodo);
               if (!ntodo) continue;
               if (td->fPartTodo->GetSize()-ntodo<500) {
                  td->fPartTodo->Set(2*td->fPartTodo->GetSize());
                  parttodo  = td->fPartTodo->GetArray();
               }
               // Do post step actions for particles suffering a given process.
               // Surviving particles are added to the next array
               propagator->Process(iproc)->PostStep(td->fVolume, ntodo, parttodo, ntotnext, partnext,tid);
               if (td->fPartNext->GetSize()-ntotnext<500) {
                  td->fPartNext->Set(2*td->fPartNext->GetSize());
                  partnext  = td->fPartNext->GetArray();
                  td->fPartInd->Set(2*td->fPartInd->GetSize());
                  particles = td->fPartInd->GetArray();
               }
            }
            memcpy(particles, partnext, ntotnext*sizeof(Int_t));
            ntotransport = ntotnext;
         }
         // I/O: Dump current generation
         //         Printf("   ### Generation %d:  %d tracks  cputime=%f", generation, ntotransport,cputime);
         if (propagator->fFillTree) {
            cputime = gPropagator->fTimer->CpuTime();
            gPropagator->fOutput->SetStamp(gPropagator->fThreadData[tid]->fVolume->GetNumber(), propagator->fWMgr->GetBasketGeneration(), generation, ntotransport, cputime);
            for (itrack=0; itrack<ntotransport;itrack++) {
               track = tracks[particles[itrack]];
               gPropagator->fOutput->SetTrack(itrack, track);
            }   
            gPropagator->fOutTree->Fill();
         }
      }   
      
   finish:
      basket->Clear();
      wm->EmptyQueue()->push(basket);
      propagator->InjectCollection(tid);
   }
   wm->AnswerQueue()->push(0);
   Printf("=== Thread %d: exiting ===", tid);
   return 0;
}
