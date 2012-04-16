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
   while ((vol=(TGeoVolume*)next())) {
      basket = new GeantVolumeBasket(vol);
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
void WorkloadManager::InterruptBasket(GeantVolumeBasket *basket, Int_t *trackin, Int_t ntracks, Int_t tid)
{
// Interrupt transporting a given basket. This has to be correlated with the fact
// that there are too many idle workers and the work queue is not filling.
   for (Int_t itr=0; itr<ntracks; itr++) basket->AddTrack(trackin[itr]);
//   basket->Clear();
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
   TThread *t = new TThread(WorkloadManager::GarbageCollect);
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
void *WorkloadManager::GarbageCollect(void *)
{
// Garbage collector thread, called by a single thread.
   GeantPropagator *propagator = GeantPropagator::Instance();
   Bool_t graphics = propagator->fUseGraphics;
   WorkloadManager *wm = WorkloadManager::Instance();
   concurrent_queue *feeder_queue = wm->FeederQueue();
   concurrent_queue *empty_queue = wm->EmptyQueue();
   // Number of baskets in the queue to transport
   Int_t ntotransport = 0;
   // Number of empty baskets available
   Int_t nempty = 0;
   Int_t ntracksperbasket = propagator->fNperBasket;
   Bool_t feed = kFALSE;
   // Feeder threshold
   Int_t min_feeder = TMath::Max(50,2*propagator->fNthreads); // setter/getter here ?
   // Number of remaining tracks
   Int_t ntracks = 1;
   // Number of tracks in the current basket
   Int_t ntracksb = 0;
   GeantBasketScheduler *basket_sch = 0;
   GeantBasketScheduler *btop = 0;
   GeantBasket *basket = 0;
   Int_t nbaskets = wm->GetNbaskets();
   Double_t factor;
   GeantVolumeBasket **array = wm->GetBasketArray();
   TH1F *hnb, *hfree;
   TCanvas *c1 = 0;
//   c1->Divide(1,2);
   TPad *pad1=0, *pad2=0;
   TFile *fresult = 0;
   if (graphics) {
      c1 = new TCanvas("c2","c2",800,900);
      c1->Divide(1,2);
      pad1 = (TPad*)c1->cd(1);
      hnb = new TH1F("hnb","number of baskets in the transport queue",500,0,500);
      hnb->SetFillColor(kRed);
      hnb->Draw();
      pad2 = (TPad*)c1->cd(2);
      hfree = new TH1F("hfree","number of free baskets available",500,0,500);
      hfree->SetFillColor(kBlue);
      hfree->Draw();
   } else {
      fresult = new TFile("results.root", "RECREATE");
      hnb = new TH1F("hnb","number of baskets in the transport queue",10000,0,10000);
      hnb->SetFillColor(kRed);
      hfree = new TH1F("hfree","average number of tracks per basket",10000,0,10000);
      hfree->SetFillColor(kBlue);
   }      
   Int_t niter = -1;
   Int_t iiter;
   Double_t nperbasket, ntracksmax;
   while (1) {
      // Monitor the queues while there are tracks to transport
      niter++;
      iiter = niter;
      nperbasket = 0;
      for (Int_t tid=0; tid<propagator->fNthreads; tid++) 
         nperbasket += propagator->fTracksPerBasket[tid];
      nperbasket /= propagator->fNthreads;
      if (graphics) iiter = niter%500;
      if (iiter==0 && graphics) {
         for (Int_t ibin=0; ibin<501; ibin++) {
            hnb->SetBinContent(ibin,0);
            hfree->SetBinContent(ibin,0);
         }  
         pad1->Modified();
         pad2->Modified();
         c1->Update();
      }   
      ntotransport = feeder_queue->size();
      nempty = empty_queue->size();
      if (iiter<10000) {
         hnb->Fill(iiter, ntotransport);
         hfree->Fill(iiter, nperbasket);
      }   
      if (graphics) {
         pad1->Modified();
         pad2->Modified();
         c1->Update();
      }
      // Try to keep 50 baskets in the queue
      factor = (Double_t)ntotransport/(2.*propagator->fNthreads);
//      factor = 1.;
      if (factor>10.) ntracksperbasket += 0.2*ntracksperbasket;
      else if (factor<0.5) ntracksperbasket -= 0.2*ntracksperbasket;
      ntracksperbasket = TMath::Max(5, ntracksperbasket);
//      propagator->fNperBasket = ntracksperbasket;
      ntracks = 0;
      // loop all baskets and garbage collect them
      if (ntotransport < min_feeder) {
//         if (!feed) Printf("=== Garbage collector: start feeding ===");
         feed=kTRUE;
      } else {
//         if (feed) Printf("=== Garbage collector: stop feeding ===");
         feed = kFALSE;
      }
      ntracksmax = 0;   
      for (Int_t ibasket=0; ibasket<nbaskets; ibasket++) {
         basket_sch = array[ibasket]->GetScheduler();
         ntracksb = basket_sch->GetNtotal();
         if (!ntracksb) continue;
         if (ntracksb>ntracksmax) {
            ntracksmax = ntracksb;
            btop = basket_sch;
         }   
         ntracks += ntracksb;
         if (!basket) basket = empty_queue->wait_and_pop();
         basket = basket_sch->GarbageCollect(basket, feed);
      }
      if ((niter%100) == 0 && btop) Printf("== Collector iteration #%d: ntracks=%d TOP=%s",
                niter, ntracks, btop->GetVolume()->GetName());
      if (!ntotransport && !ntracks) break;
   }
   Printf("=== Garbage collector: stopping threads ===");
   for (Int_t i=0; i<propagator->fNthreads; i++) wm->FeederQueue()->push(0);
   wm->AnswerQueue()->push(0);
   if (!graphics) {
      hnb->Write();
      hfree->Write();
      fresult->Close();
   }  
   Printf("=== Garbage collector: exiting ===");
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

   while (1) {
      basket = wm->FeederQueue()->wait_and_pop();
      lastToClear = kFALSE;
      if (!basket) break;
      ntotransport = basket->GetNtracks();  // all tracks to be transported 
      if (!ntotransport) goto finish;
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
      // the last thread to finish wakes up the main thread
      // ... then go to sleep
      // Checkpoint. 
//      Printf("Thread %d finished", tid);
      basket->Clear();
//      if (wm->FeederQueue()->empty()) {
//         while (!wm->IsCleared()) wm->ClearBaskets();
//         TThread::Lock();
//         while (!wm->IsFlushed()) {
//            wm->SetFlushed(wm->GetBuffer()->FlushBaskets());
//         }   
//         Printf("worker %d finished. Queue empty.", tid);
//         TThread::UnLock();
//      }
      wm->EmptyQueue()->push(basket);
//      bufferStop->StartN();
   }
   wm->AnswerQueue()->push(0);
   Printf("=== Thread %d: exiting ===", tid);
   return 0;
}
