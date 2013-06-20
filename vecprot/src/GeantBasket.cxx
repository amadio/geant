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

ClassImp(GeantBasket)

//______________________________________________________________________________
GeantBasket::GeantBasket()
            :TObject(),
             fNtracks(0),
             fMaxTracks(0),
             fIndex(0)
{
// ctor.
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t maxtracks)
            :TObject(),
             fNtracks(0),
             fMaxTracks(maxtracks),
             fIndex(new Int_t[maxtracks])
{
// ctor.
}

//______________________________________________________________________________
GeantBasket::~GeantBasket()
{
// dtor.
   delete [] fIndex;
}
   
//______________________________________________________________________________
void GeantBasket::AddTrack(Int_t itrack)
{
// Add a new track to this basket;
   if (fNtracks >= fMaxTracks-1) Resize(2*fMaxTracks);
   fIndex[fNtracks++] = itrack;
}

//______________________________________________________________________________
void GeantBasket::AddTracks(const Int_t *array, Int_t ntracks)
{
// Add array of tracks to the basket.
   if (fNtracks+ntracks > fMaxTracks-1) Resize(TMath::Max(fNtracks+ntracks, 2*fMaxTracks));
   memcpy(&fIndex[fNtracks], array, ntracks*sizeof(Int_t));
   fNtracks += ntracks;
}
   
//______________________________________________________________________________
void GeantBasket::Clear(Option_t *)
{
// Clear basket;
   fNtracks = 0;
}   

//______________________________________________________________________________
Bool_t GeantBasket::Contains(Int_t event) const
{
// Checks if any of the array of tracks belongs to the given event.
   for (Int_t itr=0; itr<fNtracks; itr++) {
      if (gPropagator->fTracks[fIndex[itr]]->event == event) return kTRUE;
   }
   return kFALSE;
}      

//______________________________________________________________________________
void GeantBasket::Print(Option_t *) const
{
// Print basket content.
   TThread::Lock();
   TString s = Form("basket %p with %d tracks:", this, fNtracks);
   for (Int_t i=0; i<fNtracks; i++) {
//      gPropagator->fTracks[fIndex[i]]->Print();
      s += Form(" %d",fIndex[i]);
   }
   Printf("%s", s.Data());
   TThread::UnLock();
}

//______________________________________________________________________________
void GeantBasket::Resize(Int_t newSize)
{
// Resize the array of track indices. Not thread safe - should be called by a 
// single thread at a time;
   Int_t *newindex = new Int_t[newSize];
   memcpy(newindex, fIndex, fNtracks*sizeof(Int_t));
   delete [] fIndex;
   fIndex = newindex;
   fMaxTracks = newSize;
}
   
ClassImp(GeantTrackCollection)

//______________________________________________________________________________
GeantTrackCollection::GeantTrackCollection()
                    :TObject(),
                     fNtracks(0),
                     fSize(0),
                     fTracks(0),
                     fBaskets(0)
{
// Default ctor.
}

//______________________________________________________________________________
GeantTrackCollection::GeantTrackCollection(Int_t size)
                    :TObject(),
                     fNtracks(0),
                     fSize(size),
                     fTracks(0),
                     fBaskets(0)
{
// Ctor.
   fTracks = new Int_t[size];
   fBaskets = new GeantVolumeBasket*[size];
}

//______________________________________________________________________________
GeantTrackCollection::~GeantTrackCollection()
{
// Dtor.
   delete [] fTracks;
   delete [] fBaskets;
}
   
//______________________________________________________________________________
Int_t GeantTrackCollection::AddTrack(Int_t itrack, GeantVolumeBasket *basket)
{
// Add a new track entering the basket.
   if (!gPropagator->fTracks[itrack]->IsAlive()) return fNtracks;
   if (fNtracks==fSize-1) {
      Int_t *tracks = new Int_t[2*fSize];
      GeantVolumeBasket **baskets = new GeantVolumeBasket*[2*fSize];
      memcpy(tracks, fTracks, fSize*sizeof(Int_t));
      memcpy(baskets, fBaskets, fSize*sizeof(GeantVolumeBasket*));
      delete [] fTracks; fTracks = tracks;
      delete [] fBaskets; fBaskets = baskets;
      fSize *= 2;
   }
   fTracks[fNtracks] = itrack;
   fBaskets[fNtracks] = basket;
   return fNtracks++;
}

//______________________________________________________________________________
Int_t GeantTrackCollection::FlushTracks(GeantMainScheduler *main)
{
// Flush all tracks to the main scheduler. Returns number of injected baskets.
   Int_t ninjected = 0;
   for (Int_t itr=0; itr<fNtracks; itr++) ninjected += main->AddTrack(fTracks[itr], fBaskets[itr]->GetNumber());
   fNtracks = 0;
   return ninjected;
}

//______________________________________________________________________________
void GeantTrackCollection::Print(Option_t *) const
{
// Print info
   TThread::Lock();
   TString s = Form("collection %p with %d tracks:", this, fNtracks);
   for (Int_t i=0; i<fNtracks; i++) {
      s += Form(" %d",fTracks[i]);
   }
   Printf("%s", s.Data());
   TThread::UnLock();
}   

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
