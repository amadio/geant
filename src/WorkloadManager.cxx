#include "WorkloadManager.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
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
#include "PhysicsProcess.h"

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
void WorkloadManager::Print()
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
   for (Int_t ith=0; ith<fNthreads; ith++) {
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
   // Flag for event dumping being effective
   Bool_t dump_effective = kFALSE;
   // Event range to be flushed
   Int_t evt_range[2] = {-1, -1};
//   GeantVolumeBasket **array = wm->GetBasketArray();
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
//   Double_t ntracksmax;
   GeantTrackCollection *collector;
   TObject **carray = new TObject*[500];
//   Bool_t direct_feed = kFALSE;
   while ((collector = (GeantTrackCollection*)collector_queue->wait_and_pop_max(500,npop,carray))) {
      // Monitor the queues while there are tracks to transport
      niter++;
      // Check first the queue of collectors
      ncollectors = npop+collector_queue->size_async();
//      Printf("Popped %d collections, %d still in the queue", npop, ncollectors-npop);
      nempty = empty_queue->size_async();
      // Process popped collectors and flush their tracks
      ninjected = 0;
      for (UInt_t icoll=0; icoll<npop; icoll++) {
         collector = (GeantTrackCollection*)carray[icoll];
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
            default:
               break;
         }
      }
      if (evt_range[0] >= 0) {
         npriority = main_sch->GetNpriority();
         if (npriority) {
            if (!dump_effective) dump_effective = kTRUE;
            ninjected += main_sch->FlushPriorityBaskets();
//            Printf("Number of priority baskets flushed: %d", npriority);
         } else {
            if (dump_effective && evt_range[0]>=0) {
               for (Int_t iev=evt_range[0]; iev<=evt_range[1]; iev++) {
                  GeantEvent *evt = propagator->fEvents[iev];
                  if (!evt->Transported()) break;
                  if (iev == evt_range[1]) {
                     dump_effective = kFALSE;
                     Printf("+++ Events %d-%d transported +++", evt_range[0],evt_range[1]);
                     evt_range[0] = evt_range[1]+1;
                     evt_range[1] = evt_range[0]+4;
                     if (evt_range[0]<propagator->fNevents-1) {
                        if (evt_range[1]>=propagator->fNevents) evt_range[1] = propagator->fNevents-1;
                        Printf("Setting priority range %d to %d", evt_range[0],evt_range[1]);
                     } else {
                        evt_range[0] = evt_range[1] = -1;
                     }   
                     main_sch->SetPriorityRange(evt_range[0], evt_range[1]);
                  }
               }
            }
         }   
      }
      ntotransport = feeder_queue->size_async();
      nwaiting = propagator->GetNwaiting();
      // If no collectors and to few baskets to transport perform a garbage collection
//      Printf("picked=%d ncoll=%d  ninjected=%d ntotransport=%d",npop, ncollectors,ninjected,ntotransport);
     if (ntotransport<min_feeder) {
        crtphase = 1;
        ninjected += main_sch->FlushBaskets();
//        Printf("Garbage collection injected %d baskets", ninjected);
        if (evt_range[0]<0) {
           // Start dumping
           evt_range[0] = 0;
           evt_range[1] = 4;
           Printf("Setting priority range %d to %d", evt_range[0],evt_range[1]);
           main_sch->SetPriorityRange(evt_range[0], evt_range[1]);
        }
        nwaiting = propagator->GetNwaiting();
        if (!ninjected && !ntotransport && nwaiting==nworkers) break;
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
   for (Int_t i=0; i<propagator->fNevents; i++) {
      GeantEvent *evt = propagator->fEvents[i];
      Printf("Event %d: total=%d transported=%d", i, evt->ntracks, evt->ndone);
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
   for (Int_t i=0; i<propagator->fNthreads; i++) wm->FeederQueue()->push(0);
   wm->AnswerQueue()->push(0);
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
   Int_t ntotnext, ntmp, ntodo, ncross, cputime, ntotransport;
   GeantTrack *track = 0;
   Int_t itrack;
   Int_t generation = 0;
   Bool_t lastToClear = kFALSE;
   GeantVolumeBasket *basket_sch = 0;
   GeantBasket *basket = 0;
   Int_t tid = TGeoManager::ThreadId();
   GeantPropagator *propagator = GeantPropagator::Instance();
   WorkloadManager *wm = WorkloadManager::Instance();
   Int_t *particles = propagator->fPartInd[tid]->GetArray();
   Int_t *partnext  = propagator->fPartNext[tid]->GetArray();
   Int_t *parttodo  = gPropagator->fPartTodo[tid]->GetArray();
   Int_t *partcross = gPropagator->fPartCross[tid]->GetArray();
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
      propagator->fWaiting[tid] = 1;
      basket = (GeantBasket*)wm->FeederQueue()->wait_and_pop();
      propagator->fWaiting[tid] = 0;
      lastToClear = kFALSE;
      if (!basket) break;
      ntotransport = basket->GetNtracks();  // all tracks to be transported 
      if (!ntotransport) goto finish;
//      Printf("======= BASKET taken by thread #%d =======", tid);
//      basket->Print();
//      Printf("==========================================");
      propagator->fTracksPerBasket[tid] = ntotransport;
      path = tracks[basket->GetTracks()[0]]->path;
      gPropagator->fVolume[tid] = path->GetCurrentNode()->GetVolume();
      basket_sch = (GeantVolumeBasket*)gPropagator->fVolume[tid]->GetField();
      wm->SetCurrentBasket(tid,basket_sch);
//      Printf("(%d) ================= BASKET of %s: %d tracks", tid, gPropagator->fVolume[tid]->GetName(), ntotransport);
      memcpy(particles, basket->GetTracks(), ntotransport*sizeof(Int_t));
//      PrintParticles(particles, ntotransport, tid);
//      sprintf(slist,"Thread #%d transporting %d tracks in basket %s: ", tid, ntotransport, basket->GetName());
//      sslist = slist;
//      for (Int_t ip=0; ip<ntotransport; ip++) {sprintf(slist,"%d ",particles[ip]), sslist += slist;}
//      Printf("%s", sslist.Data());
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
//            printf("(%d) %s   ->crossing particles (%d): ",tid, basket->GetName(), ncross); 
//            for (Int_t ii=0; ii<ncross; ii++) printf("%d ", partcross[ii]);
//            printf("\n(%d) %s   ->remaining particles (%d): ", tid, basket->GetName(), ntotnext);
//            for (Int_t ii=0; ii<ntotnext; ii++) printf("%d ", partnext[ii]);
//            printf("\n(%d) %s   ->todo particles (%d): ", tid, basket->GetName(),ntodo);
//            for (Int_t ii=0; ii<ntodo; ii++) printf("%d ", parttodo[ii]);
//            printf("\n");
            // Post-step actions by continuous processes for particles reaching boundaries
            if (propagator->fUsePhysics && ncross) {
               for (Int_t iproc=0; iproc<nprocesses; iproc++) {
                  if (propagator->Process(iproc)->IsType(PhysicsProcess::kDiscrete)) continue;
                  Int_t nafter = 0;
                  gPropagator->Process(iproc)->PostStep(gPropagator->fVolume[tid], ncross, partcross, nafter, NULL,tid);
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
               gPropagator->Process(iproc)->PostStep(gPropagator->fVolume[tid], ntotransport, particles, ntodo, parttodo, tid);
               // Do we have stopped particles ?
               if (ntodo<ntotransport) {
                  memcpy(particles, parttodo, ntodo*sizeof(Int_t));
                  ntotransport = ntodo;
               }
            } 
            // Copy al tracks for which step was limited by a continuous process
            // to the next array
            for (Int_t itr=0; itr<ntotransport; itr++) {
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
               if (gPropagator->fPartTodo[tid]->GetSize()-ntodo<500) {
                  gPropagator->fPartTodo[tid]->Set(2*gPropagator->fPartTodo[tid]->GetSize());
                  parttodo  = gPropagator->fPartTodo[tid]->GetArray();
               }   
               // Do post step actions for particles suffering a given process.
               // Surviving particles are added to the next array
               gPropagator->Process(iproc)->PostStep(gPropagator->fVolume[tid], ntodo, parttodo, ntotnext, partnext,tid);
               if (gPropagator->fPartNext[tid]->GetSize()-ntotnext<500) {
                  gPropagator->fPartNext[tid]->Set(2*gPropagator->fPartNext[tid]->GetSize());
                  partnext  = gPropagator->fPartNext[tid]->GetArray();
                  propagator->fPartInd[tid]->Set(2*propagator->fPartInd[tid]->GetSize());
                  particles = propagator->fPartInd[tid]->GetArray();
               }   
            }
            memcpy(particles, partnext, ntotnext*sizeof(Int_t));
            ntotransport = ntotnext;
         }
         // I/O: Dump current generation
//         Printf("   ### Generation %d:  %d tracks  cputime=%f", generation, ntotransport,cputime);
         if (propagator->fFillTree) {
            cputime = gPropagator->fTimer->CpuTime();
            gPropagator->fOutput->SetStamp(gPropagator->fVolume[tid]->GetNumber(), propagator->fWMgr->GetBasketGeneration(), generation, ntotransport, cputime);
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
