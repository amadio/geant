#include "GeantBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "PhysicsProcess.h"
#include "WorkloadManager.h"

#include "TThread.h"
#include "TArrayI.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TGeoBranchArray.h"

ClassImp(GeantMainScheduler)

//______________________________________________________________________________
GeantMainScheduler::GeantMainScheduler()
                   :TObject(),
                    fNvolumes(0),
                    fNpriority(0),
                    fBaskets(0),
                    fPriorityBaskets(0),
                    feeder_queue(0),
                    empty_queue(0),
                    collector_queue(0)
{
// dummy
   fPriorityRange[0] = fPriorityRange[1] = -1;
}

//______________________________________________________________________________
GeantMainScheduler::GeantMainScheduler(Int_t nvolumes)
                   :TObject(),
                    fNvolumes(nvolumes),
                    fNpriority(0),
                    fBaskets(new GeantBasket*[nvolumes]),
                    fPriorityBaskets(new GeantBasket*[nvolumes]),
                    feeder_queue(0),
                    empty_queue(0),
                    collector_queue(0)
{
// ctor
   WorkloadManager *wm = WorkloadManager::Instance();
   feeder_queue = wm->FeederQueue();
   empty_queue = wm->EmptyQueue();
   collector_queue = wm->CollectorQueue();
   for (Int_t i=0; i<nvolumes; i++) {
      fBaskets[i] = new GeantBasket(2*gPropagator->fNperBasket);
      fPriorityBaskets[i] = new GeantBasket(2*gPropagator->fNperBasket);
   }         
   fPriorityRange[0] = fPriorityRange[1] = -1;
}

//______________________________________________________________________________
GeantMainScheduler::~GeantMainScheduler()
{
// dtor.
   if (fBaskets) {
      for (Int_t ib=0; ib<fNvolumes; ib++) {
         delete fBaskets[ib];
         delete fPriorityBaskets[ib];
      }   
   }   
   delete [] fBaskets;                    
   delete [] fPriorityBaskets;                    
}

//______________________________________________________________________________
Int_t GeantMainScheduler::AddTrack(GeantTrack_v &tracks)
{
// Add all tracks and inject baskets if above threshold. Returns the number
// of injected baskets. Baskets pushed at fNperBasket.
   Int_t ninjected = 0;
   Bool_t priority = kFALSE;
   GeantBasket **baskets = fBaskets;
   Int_t ntracks = tracks.GetNtracks();
   for (Int_t itr=0; itr<ntracks; itr++) {
      if (tracks.fStatusV[itr]==GeantTrack::kKilled) continue;
      if (fPriorityRange[0]>=0) {
         Int_t event = tracks.fEventV[itr];
         priority = kFALSE;
         if (event>=fPriorityRange[0] && event<=fPriorityRange[1]) {
            baskets = fPriorityBaskets;
            priority = kTRUE;
            if (!baskets[ibasket]->GetNtracks()) fNpriority++;
         }
      }
      Int_t level = tracks.fPathV[itr]->GetLevel();
      Int_t ibasket = tracks.fPathV[itr]->GetNode(level)->GetVolume()->GetNumber();
      baskets[ibasket]->AddTrack(tracks, itr);
      if (baskets[ibasket]->GetNtracks() >= gPropagator->fNperBasket) {
   // inject this basket
//      Printf("pushing basket %p", baskets[ibasket]);
//      baskets[ibasket]->Print();
         feeder_queue->push(baskets[ibasket], priority);
         fNpriority -= (Int_t)priority;
         ninjected++;
         if (empty_queue->empty_async()) baskets[ibasket] = new GeantBasket(2*gPropagator->fNperBasket);
         else {
            baskets[ibasket] = (GeantBasket*)empty_queue->wait_and_pop();
            baskets[ibasket]->Clear();
         }
      }
   }   
   return ninjected;
}  
   
//______________________________________________________________________________
Int_t GeantMainScheduler::AddTrack(Int_t itrack, Int_t ibasket)
{
// Add track and inject basket if above threshold. Returns 1 if the basket was
// injected and 0 if not. Basket pushed at fNperBasket
   Int_t ninjected = 0;
   Bool_t priority = kFALSE;
   GeantBasket **baskets = fBaskets;
   GeantTrack *track = gPropagator->fTracks[itrack];
   if (!track || !track->IsAlive()) {
//      Printf("TRACK %d is NOT ALIVE", itrack);
      return 0;
   }  
   if (fPriorityRange[0]>=0) {
      Int_t event = gPropagator->fTracks[itrack]->event;
      priority = kFALSE;
      if (event>=fPriorityRange[0] && event<=fPriorityRange[1]) {
         baskets = fPriorityBaskets;
         priority = kTRUE;
         if (!baskets[ibasket]->GetNtracks()) fNpriority++;
      }   
   }
   baskets[ibasket]->AddTrack(itrack);
   if (baskets[ibasket]->GetNtracks() >= gPropagator->fNperBasket) {
   // inject this basket
//      Printf("pushing basket %p", baskets[ibasket]);
//      baskets[ibasket]->Print();
      feeder_queue->push(baskets[ibasket], priority);
      fNpriority -= (Int_t)priority;
      ninjected++;
      if (empty_queue->empty_async()) baskets[ibasket] = new GeantBasket(2*gPropagator->fNperBasket);
      else {
         baskets[ibasket] = (GeantBasket*)empty_queue->wait_and_pop();
         baskets[ibasket]->Clear();
      }
   }
   return ninjected;
}

//______________________________________________________________________________
Int_t GeantMainScheduler::FlushPriorityBaskets()
{
// Flush all non-empty priority baskets.
   if (!fNpriority) return 0;
   Int_t ninjected = 0;
   for (Int_t ibasket=0; ibasket<fNvolumes; ibasket++) {
      if (!fPriorityBaskets[ibasket]->GetNtracks()) continue;
      // inject this basket
//      Printf("pushing basket %p", fPriorityBaskets[ibasket]);
//      fPriorityBaskets[ibasket]->Print();
      feeder_queue->push(fPriorityBaskets[ibasket], kTRUE);
      fNpriority--;
      ninjected++;
      if (empty_queue->empty_async()) fPriorityBaskets[ibasket] = new GeantBasket(2*gPropagator->fNperBasket);
      else {
         fPriorityBaskets[ibasket] = (GeantBasket*)empty_queue->wait_and_pop();
         fPriorityBaskets[ibasket]->Clear();
      }   
   }
   return ninjected;
}   

//______________________________________________________________________________
Int_t GeantMainScheduler::FlushBaskets(Int_t threshold)
{
// Flush baskets in the work queue until the number of injected objects is above
// threshold. If event is specified, flush only baskets containing the event.
   Int_t ntoflush = fNvolumes;
   Int_t ninjected = 0;
   if (threshold) ntoflush = threshold;
   Int_t ntotransport = feeder_queue->size_async();
   ntoflush -= ntotransport;
   if (ntoflush < 0) return 0;
   for (Int_t ibasket=0; ibasket<fNvolumes; ibasket++) {
      // inject this basket
      if (!fBaskets[ibasket]->GetNtracks()) continue;
//      Printf("pushing basket %p", fBaskets[ibasket]);
//      fBaskets[ibasket]->Print();
      feeder_queue->push(fBaskets[ibasket], kFALSE);
      ninjected++;
      ntoflush--;
      if (empty_queue->empty()) fBaskets[ibasket] = new GeantBasket(2*gPropagator->fNperBasket);
      else {
         fBaskets[ibasket] = (GeantBasket*)empty_queue->wait_and_pop();
         fBaskets[ibasket]->Clear();
      }   
      if (ntoflush<=0) return ninjected;
   }
   return ninjected;
}
