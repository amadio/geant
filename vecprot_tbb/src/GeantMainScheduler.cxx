#include "GeantMainScheduler.h"

#include "GeantPropagator.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "GeantTrack.h"

ClassImp(GeantMainScheduler)

GeantMainScheduler* GeantMainScheduler::fgInstance = 0;

//______________________________________________________________________________
GeantMainScheduler::GeantMainScheduler()
                   :TObject(),
                    fNvolumes(0),
                    fNpriority(0),
                    fBaskets(0),
                    fPriorityBaskets(0)
{
// dummy
   fgInstance = this;
}

//______________________________________________________________________________
GeantMainScheduler::GeantMainScheduler(Int_t nvolumes)
                   :TObject(),
                    fNvolumes(nvolumes),
                    fNpriority(0),
                    fBaskets(new GeantBasket*[nvolumes]),
                    fPriorityBaskets(new GeantBasket*[nvolumes])
{
// ctor
   GeantPropagator* propagator = GeantPropagator::Instance();
   for (Int_t i=0; i<nvolumes; i++) {
      fBaskets[i] = new GeantBasket(2*propagator->fNperBasket);
      fPriorityBaskets[i] = new GeantBasket(2*propagator->fNperBasket);
   }
   fgInstance = this;
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

GeantMainScheduler* GeantMainScheduler::Instance(Int_t nvolumes)
{
// Return singleton instance.
   if (fgInstance) return fgInstance;
   if (!nvolumes) {
      ::Error("GeantMainScheduler::Instance", "No instance yet so you should provide number of volumes.");
      return 0;
   }
   return new GeantMainScheduler(nvolumes);
}

//______________________________________________________________________________
Int_t GeantMainScheduler::AddTrack(Int_t itrack, Int_t ibasket, Bool_t* pushedPriority)
{
// Add track and inject basket if above threshold. Returns 1 if the basket was
// injected and 0 if not. Basket pushed at fNperBasket

   GeantPropagator *gPropagator = GeantPropagator::Instance();
   WorkloadManager *wm = WorkloadManager::Instance();

   Int_t ninjected = 0;
   Bool_t priority = kFALSE;
   GeantBasket **baskets = fBaskets;
   GeantTrack *track = gPropagator->fTracks[itrack];
   if (!track || !track->IsAlive()) {
//      Printf("TRACK %d is NOT ALIVE", itrack);
      return 0;
   }
   if (gPropagator->fPrioritize) {
      Int_t event = gPropagator->fTracks[itrack]->event;
      if (event>=gPropagator->fPriorityRange[0] && event<=gPropagator->fPriorityRange[1]) {
         baskets = fPriorityBaskets;
         priority = kTRUE;
         if (!baskets[ibasket]->GetNtracks()) fNpriority++;
      }
   }

   // CRITICAL SECTION BEGIN
   the_lock.Lock();

   baskets[ibasket]->AddTrack(itrack);
   if (baskets[ibasket]->GetNtracks() >= gPropagator->fNperBasket 
         /*&& wm->feeder_queue->size_async() < 1000 ||
         baskets[ibasket]->GetNtracks() <= gPropagator->fMaxPerBasket/2*/) {
   // inject this basket
//      Printf("pushing basket %p", baskets[ibasket]);
//      baskets[ibasket]->Print();


      if (!priority) wm->tbb_feeder_queue.push(baskets[ibasket]);
      else wm->tbb_feeder_priority_queue.push(baskets[ibasket]);

      fNpriority -= (Int_t)priority;
      ninjected++;
      *pushedPriority = priority;

      if (wm->tbb_feeder_empty_queue.try_pop(baskets[ibasket])) baskets[ibasket]->Clear();
      else baskets[ibasket] = new GeantBasket(2*gPropagator->fNperBasket);

   }

   the_lock.UnLock();
   // CRITICAL SECTION END

   return ninjected;
}

//______________________________________________________________________________
Int_t GeantMainScheduler::FlushPriorityBaskets()
{
// Flush all non-empty priority baskets.

   GeantPropagator *gPropagator = GeantPropagator::Instance();
   WorkloadManager *wm = WorkloadManager::Instance();

   if (!fNpriority) return 0;
   Int_t ninjected = 0;

   // CRITICAL SECTION BEGIN
   the_lock.Lock();

   for (Int_t ibasket=0; ibasket<fNvolumes; ibasket++) {
      if (!fPriorityBaskets[ibasket]->GetNtracks()) continue;
      // inject this basket
//      Printf("pushing basket %p", fPriorityBaskets[ibasket]);
//      fPriorityBaskets[ibasket]->Print();

      wm->tbb_feeder_priority_queue.push(fPriorityBaskets[ibasket]);

      fNpriority--;
      ninjected++;

      if (wm->tbb_feeder_empty_queue.try_pop(fPriorityBaskets[ibasket])) fPriorityBaskets[ibasket]->Clear();
      else fPriorityBaskets[ibasket] = new GeantBasket(2*gPropagator->fNperBasket);

   }

   the_lock.UnLock();
   // CRITICAL SECTION END

   return ninjected;
}

//______________________________________________________________________________
Int_t GeantMainScheduler::FlushNormalBaskets(Int_t threshold)
{
// Flush baskets in the work queue until the number of injected objects is above
// threshold. If event is specified, flush only baskets containing the event.

   GeantPropagator *gPropagator = GeantPropagator::Instance();
   WorkloadManager *wm = WorkloadManager::Instance();

   Int_t ntoflush = fNvolumes;
   Int_t ninjected = 0;
   if (threshold) ntoflush = threshold;
   Int_t ntotransport = wm->tbb_feeder_queue.size();

   ntoflush -= ntotransport;
   if (ntoflush < 0) return 0;

   // CRITICAL SECTION BEGIN
   the_lock.Lock();

   for (Int_t ibasket=0; ibasket<fNvolumes; ibasket++) {
      // inject this basket
      if (!fBaskets[ibasket]->GetNtracks()) continue;
//      Printf("pushing basket %p", fBaskets[ibasket]);
//      fBaskets[ibasket]->Print();

      wm->tbb_feeder_queue.push(fBaskets[ibasket]);

      ninjected++;
      ntoflush--;
      if (wm->tbb_feeder_empty_queue.try_pop(fBaskets[ibasket])) fBaskets[ibasket]->Clear();
      else fBaskets[ibasket] = new GeantBasket(2*gPropagator->fNperBasket);


      if (ntoflush<=0) return ninjected;
   }

   the_lock.UnLock();
   // CRITICAL SECTION END

   return ninjected;
}

