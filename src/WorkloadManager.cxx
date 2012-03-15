#include "WorkloadManager.h"

#include "TList.h"
#include "TTree.h"
#include "TThread.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "GeantVolumeBasket.h"
#include "GeantParticleBuffer.h"
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
   fNminThreshold = 50;
   fNqueued = 0;
   fBindex = 0;
   fStarted = kFALSE;
   feeder_queue = new concurrent_queue(true);
   answer_queue = new concurrent_queue;
   fgInstance = this;
   fCurrentBasket = new GeantVolumeBasket*[nthreads];
   fListThreads = 0;
   fBasketArray = 0;
   fBuffer = new GeantParticleBuffer(nthreads, 100000);
}

//______________________________________________________________________________
WorkloadManager::~WorkloadManager()
{
// Destructor.
   delete feeder_queue;
   delete answer_queue;
   delete [] fCurrentBasket;
   if (fBasketArray) {
      for (Int_t i=0; i<fNbaskets; i++) delete fBasketArray[i];
      delete [] fBasketArray;
   }   
   delete [] fBindex;
   delete fBuffer;
   fgInstance = 0;
}

//______________________________________________________________________________
void WorkloadManager::AddPendingTrack(Int_t itrack, GeantVolumeBasket *basket, Int_t tid)
{
// Add a pending track.
   fBuffer->AddPendingTrack(itrack, basket, tid);
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
   fBindex = new Int_t[nvolumes];
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
void WorkloadManager::ClearBaskets()
{
// Clear all active particles from the transported baskets. Also flush all 
// pending particles
   for (Int_t ibasket=0; ibasket<fNbasketgen; ibasket++)
      fBasketArray[fBindex[ibasket]]->Clear();
   fBuffer->FlushBaskets();   
}   

//______________________________________________________________________________
void WorkloadManager::Print()
{
//
   feeder_queue->Print();
}   

//______________________________________________________________________________
void WorkloadManager::QueueBaskets()
{
// Queue all baskets in the current generation. Note that workers will start as
// soon as the first chunk gets pushed in the queue.
   GeantVolumeBasket *basket = 0;
   Int_t nchunks;
   for (Int_t ibasket=0; ibasket<fNbasketgen; ibasket++) {
      basket = fBasketArray[fBindex[ibasket]];
//      basket->Prepare();
      nchunks = basket->GetNchunks(fNthreads);
      for (Int_t ichunk=0; ichunk<nchunks; ichunk++) {
         fNqueued++;
         feeder_queue->push(basket);
      }   
   }
}   
   
//______________________________________________________________________________
void WorkloadManager::SelectBaskets()
{
// Select the list of baskets to be transported in the current generation.
// Called by main thread.
//   SortBaskets();
//   Bool_t useThreshold = kFALSE;
//   if (fBasketArray[fBindex[0]]->GetNtotal()>fNminThreshold) useThreshold = kTRUE;
   fNbasketgen = 0;
   Int_t ntrackgen = 0;
   Int_t indmax = 0;
   Int_t nmax = 0;
   for (Int_t ibasket=0; ibasket<fNbaskets; ibasket++) {
      Int_t ntracks = fBasketArray[ibasket]->GetNtotal();
      if (ntracks<fNminThreshold) continue;
      fBindex[fNbasketgen++] = ibasket;
      ntrackgen += ntracks;
      if (ntracks>nmax) {
         nmax = ntracks;
         indmax = ibasket;
      }   
   }
   if (!fNbasketgen) {
      Printf("Garbage collection");
      for (Int_t ibasket=0; ibasket<fNbaskets; ibasket++) {
         Int_t ntracks = fBasketArray[ibasket]->GetNtotal();
         if (!ntracks) continue;
         if (ntracks>nmax) {
            nmax = ntracks;
            indmax = ibasket;
         }
         fBindex[fNbasketgen++] = ibasket;
         ntrackgen += ntracks;
      }   
   }
   if (!fStarted) StartThreads();
   if (fNbaskets) Printf("#### GENERATION #%04d (%05d part) OF %04d/%04d VOLUME BASKETS, (TOP= %s - %d tracks) ####", fBasketGeneration, ntrackgen, fNbasketgen, fNbaskets, fBasketArray[indmax]->GetName(), nmax);
}   

//______________________________________________________________________________
void WorkloadManager::SortBaskets()
{
// Sort baskets in decreasing number of particles. The order is set in the provided index array of size fNbaskets minimum.
   if (fNbaskets<1) return;
   Int_t *array = new Int_t[fNbaskets];
   array[0] = 0;
   for (Int_t ibasket=0; ibasket<fNbaskets; ibasket++) array[ibasket] = fBasketArray[ibasket]->GetNtotal();
   TMath::Sort(fNbaskets, array, fBindex, kTRUE);
   delete [] array;
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
}   

//______________________________________________________________________________
void WorkloadManager::JoinThreads()
{
// 
   for (Int_t ith=0; ith<fNthreads; ith++) feeder_queue->push(0);
   for (Int_t ith=0; ith<fNthreads; ith++) ((TThread*)fListThreads->At(ith))->Join();
}
   
//______________________________________________________________________________
void WorkloadManager::WaitWorkers()
{
// Waiting point for the main thread until work gets done.
   while(fNqueued) {
      answer_queue->wait_and_pop();
      fNqueued--;
//      Printf("Worker %d finished", ipop);
   }
   fBasketGeneration++;
}

//______________________________________________________________________________
void *WorkloadManager::TransportTracks(void *)
{
// Thread propagating all tracks from a basket.
//      char slist[256];
//      TString sslist;
   Int_t indmin, indmax;
   Int_t ntotnext, ntmp, ntodo, ncross, cputime, ntotransport;
   GeantTrack *track = 0;
   Int_t itrack;
   Int_t *particles = 0;
   Int_t *partnext  = 0;
   Int_t *parttodo  = 0;
   Int_t *partcross = 0;
   Int_t generation = 0;
   GeantVolumeBasket *basket = 0;
   Int_t tid = TGeoManager::ThreadId();
   GeantPropagator *propagator = GeantPropagator::Instance();
   WorkloadManager *wm = WorkloadManager::Instance();
   Int_t nprocesses = propagator->fNprocesses;
   Bool_t useDebug = propagator->fUseDebug;
   GeantTrack **tracks = propagator->fTracks;
//   Printf("(%d) WORKER started", tid);
   // Create navigator if none serving this thread.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) nav = gGeoManager->AddNavigator();

   while (1) {
      basket = wm->FeederQueue()->wait_and_pop();
      if (!basket) return 0;
      wm->SetCurrentBasket(tid,basket);
//      Printf("Popped %d\n", ntmp);
//      bufferStart->Receive();
      gPropagator->fVolume[tid] = basket->GetVolume();
      ntotransport = basket->GetNtracks();  // all tracks to be transported      
      if (!ntotransport) goto finish;
      // Work splitting per thread
      basket->GetWorkload(indmin, indmax);
      ntotransport = indmax-indmin;
      if (!ntotransport) goto finish;      
//      Printf("(%d) ================= BASKET %s: %d tracks (%d-%d)", tid, basket->GetName(), ntotransport, indmin,indmax);
      particles = propagator->fPartInd[tid]->GetArray();
      partnext  = propagator->fPartNext[tid]->GetArray();
      parttodo  = gPropagator->fPartTodo[tid]->GetArray();
      partcross = gPropagator->fPartCross[tid]->GetArray();
      memcpy(particles, &basket->GetIndArray()[indmin], ntotransport*sizeof(Int_t));
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
         generation++;
         // Loop all tracks to generate physics/geometry steps
         // Physics step
         if (propagator->fUsePhysics) propagator->PhysicsSelect(ntotransport, particles, tid);
         // Geometry snext and safety
         basket->ComputeTransportLength(ntotransport, particles);
         // Propagate tracks with physics step and check if they survive boundaries 
         // or physics
         ntmp = ntotransport;
         ntodo = ntotransport;
         ntotnext = 0;
         Int_t *ptrParticles = particles;
         // Propagate all tracks alive with physics step.
         // If a boundary is encountered before the physics step, stop the track
         if (useDebug) Printf("(%d) --- propagating %d tracks for volume basket %s", tid, ntodo, basket->GetName());
         while (ntodo) {
            ntodo = 0;
            ncross = 0;
            // Propagate ALL ntmp tracks
            basket->PropagateTracks(ntmp, ptrParticles, ntotnext, partnext, ntodo, parttodo, ncross, partcross);
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
                  basket->ResetStep(ncross, partcross);
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
      //         Printf("PostStep for proc %d: %d particles:\n", iproc, ntodo);
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
            gPropagator->fOutput->SetStamp(basket->GetVolume()->GetNumber(), propagator->fWMgr->GetBasketGeneration(), generation, ntotransport, cputime);
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
      wm->AnswerQueue()->push(basket);
//      bufferStop->StartN();
   }
   return 0;
}
