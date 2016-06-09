#include "GeantMainScheduler.h"

#include "GeantPropagator.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "GeantTrack.h"

GeantMainScheduler* GeantMainScheduler::fgInstance = 0;

//______________________________________________________________________________
GeantMainScheduler::GeantMainScheduler()
                   :TObject(),
                    fNvolumes(0),
//                    fNpriority(0),
                    fBaskets(0),
                    fPriorityBaskets(0)
{
// dummy
   fgInstance = this;
//   Printf ("GeantMainScheduler dummy constructed");
}

//______________________________________________________________________________
GeantMainScheduler::GeantMainScheduler(int nvolumes)
                   :TObject(),
                    fNvolumes(nvolumes),
//                    fNpriority(0),
                    fBaskets(new GeantBasket*[nvolumes]),
                    fPriorityBaskets(new GeantBasket*[nvolumes])
{
// ctor
   GeantPropagator* propagator = GeantPropagator::Instance();
   for (int i=0; i<nvolumes; i++) {
      fBaskets[i] = new GeantBasket(2*propagator->fNperBasket);
      fPriorityBaskets[i] = new GeantBasket(2*propagator->fNperBasket);
   }
   fgInstance = this;

//   Printf ("GeantMainScheduler constructed");
}

//______________________________________________________________________________
GeantMainScheduler::~GeantMainScheduler()
{
// dtor.
   if (fBaskets) {
      for (int ib=0; ib<fNvolumes; ib++) {
         delete fBaskets[ib];
         delete fPriorityBaskets[ib];
      }
   }
   delete [] fBaskets;
   delete [] fPriorityBaskets;
}

//______________________________________________________________________________
GeantMainScheduler* GeantMainScheduler::Instance(int nvolumes)
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
int GeantMainScheduler::AddTrack(int itrack, int ibasket, bool* pushedPriority)
{
   GeantPropagator *gPropagator = GeantPropagator::Instance();
   WorkloadManager *wm = WorkloadManager::Instance();

   int ninjected = 0;
   bool priority = kFALSE;
   GeantBasket **baskets = fBaskets;
   GeantTrack *track = gPropagator->fTracks[itrack];

   if (!track || !track->IsAlive()) {
      //Printf("[GeantMainScheduler::AddTrack] TRACK %d is NOT ALIVE", itrack);
      return 0;
   }

   if (gPropagator->fPriorityRange[0] > -1) {
      int event = gPropagator->fTracks[itrack]->event;
      if (event >= gPropagator->fPriorityRange[0] && event <= gPropagator->fPriorityRange[1]) {
         baskets = fPriorityBaskets;
         priority = kTRUE;
      }
   }

   // CRITICAL SECTION BEGIN
   the_lock.Lock();

   baskets[ibasket]->AddTrack(itrack);

   if (baskets[ibasket]->GetNtracks() >= gPropagator->fNperBasket)
   {
      if (!priority) wm->tbb_feeder_queue.push(baskets[ibasket]);
      else wm->tbb_feeder_priority_queue.push(baskets[ibasket]);

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
int GeantMainScheduler::FlushPriorityBaskets()
{
   GeantPropagator *propagator = GeantPropagator::Instance();
   WorkloadManager *wm = WorkloadManager::Instance();

   int ninjected = 0;

   // CRITICAL SECTION BEGIN
   the_lock.Lock();

   for (int ibasket=0; ibasket<fNvolumes; ibasket++)
   {
      if (!fPriorityBaskets[ibasket]->GetNtracks()) continue;

      wm->tbb_feeder_priority_queue.push(fPriorityBaskets[ibasket]);
      ninjected++;

      if (wm->tbb_feeder_empty_queue.try_pop(fPriorityBaskets[ibasket])) fPriorityBaskets[ibasket]->Clear();
      else fPriorityBaskets[ibasket] = new GeantBasket(2*propagator->fNperBasket);
   }

   the_lock.UnLock();
   // CRITICAL SECTION END

   return ninjected;
}

//______________________________________________________________________________
int GeantMainScheduler::FlushNormalBaskets()
{
   GeantPropagator *propagator = GeantPropagator::Instance();
   WorkloadManager *wm = WorkloadManager::Instance();

   int ninjected = 0;

   // CRITICAL SECTION BEGIN
   the_lock.Lock();

   for (int ibasket=0; ibasket<fNvolumes; ibasket++)
   {
      if (!fBaskets[ibasket]->GetNtracks()) continue;

      wm->tbb_feeder_queue.push(fBaskets[ibasket]);
      ninjected++;

      if (wm->tbb_feeder_empty_queue.try_pop(fBaskets[ibasket])) fBaskets[ibasket]->Clear();
      else fBaskets[ibasket] = new GeantBasket(2*propagator->fNperBasket);
   }

   the_lock.UnLock();
   // CRITICAL SECTION END

   return ninjected;
}

