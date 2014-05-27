#include "TThread.h"
#include "GeantBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "GeantScheduler.h"
#include "PhysicsProcess.h"
#include "WorkloadManager.h"

#include "TThread.h"
#include "TArrayI.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNavigator.h"

ClassImp(GeantBasket)

//______________________________________________________________________________
GeantBasket::GeantBasket()
            :TObject(),
             fManager(0),
             fTracksIn(),
             fTracksOut()
{
// dummy ctor.
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, GeantBasketMgr *mgr)
            :TObject(),
             fManager(mgr),
             fTracksIn(size),
             fTracksOut(size)
{
// ctor.
}

//______________________________________________________________________________
GeantBasket::~GeantBasket()
{
// dtor.
}
   
//______________________________________________________________________________
void GeantBasket::AddTrack(GeantTrack &track)
{
// Add a new track to this basket;
   fTracksIn.AddTrack(track);
}

//______________________________________________________________________________
void GeantBasket::AddTrack(GeantTrack_v &tracks, Int_t itr)
{
// Add track from a track_v array
   fTracksIn.AddTrack(tracks, itr);
}

//______________________________________________________________________________
void GeantBasket::AddTracks(GeantTrack_v &tracks, Int_t istart, Int_t iend)
{
// Add multiple tracks from a track_v array
   fTracksIn.AddTracks(tracks, istart, iend);
}
   
//______________________________________________________________________________
void GeantBasket::Clear(Option_t *option)
{
// Clear basket;
   SetMixed(kFALSE);
   fTracksIn.Clear(option);
   fTracksOut.Clear(option);
}   

//______________________________________________________________________________
Bool_t GeantBasket::Contains(Int_t evstart, Int_t nevents) const
{
// Checks if any of the array of tracks belongs to the given event.
   return fTracksIn.Contains(evstart, nevents);
}      

//______________________________________________________________________________
TGeoVolume *GeantBasket::GetVolume() const
{
// Returns volume for this basket
   return fManager->GetVolume();
}   

//______________________________________________________________________________
void GeantBasket::Print(Option_t *) const
{
// Print basket content.
}

//______________________________________________________________________________
void GeantBasket::PrintTrack(Int_t /*itr*/, Bool_t /*input*/) const
{
// Print a given track.
}

//______________________________________________________________________________
void GeantBasket::Recycle()
{
// Recycle the basket to the volume scheduler.
   fManager->RecycleBasket(this);   
}

//______________________________________________________________________________
// Basket manager for a given volume. Holds a list of free baskets stored in a
// concurrent queue
//______________________________________________________________________________

ClassImp(GeantBasketMgr)

//______________________________________________________________________________
GeantBasketMgr::GeantBasketMgr(GeantScheduler *sch, TGeoVolume *vol, Int_t number)
                  :TGeoExtension(),
                   fScheduler(sch),
                   fVolume(vol),
                   fNumber(number),
                   fThreshold(0),
                   fNbaskets(0),
                   fNused(0),
                   fCBasket(0),
                   fPBasket(0),
                   fBaskets(),
                   fFeeder(0),
                   fMutex()
{
// Constructor
   fCBasket = GetNextBasket();
   fPBasket = GetNextBasket();
}

//______________________________________________________________________________
GeantBasketMgr::~GeantBasketMgr()
{
// Clean up
   delete fCBasket;
   delete fPBasket;
}   

//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrack(GeantTrack_v &trackv, Int_t itr, Bool_t priority)
{
// Copy directly from a track_v a track to the basket manager.
#ifdef __STAT_DEBUG
      fScheduler->GetPendingStat().AddTrack(trackv, itr);
#endif   
   if (priority) {
      fPBasket->AddTrack(trackv, itr);
      if (fPBasket->GetNinput() >= fThreshold) {
#ifdef __STAT_DEBUG
         fScheduler->GetPendingStat().RemoveTracks(fPBasket->GetInputTracks());
         fScheduler->GetQueuedStat().AddTracks(fPBasket->GetInputTracks());
#endif   
         fFeeder->push(fPBasket, priority);
         fPBasket = GetNextBasket();
         return 1;
      }
   } else {
      fCBasket->AddTrack(trackv, itr);
      if (fCBasket->GetNinput() >= fThreshold) {
#ifdef __STAT_DEBUG
         fScheduler->GetPendingStat().RemoveTracks(fCBasket->GetInputTracks());
         fScheduler->GetQueuedStat().AddTracks(fCBasket->GetInputTracks());
#endif   
         fFeeder->push(fCBasket, priority);
         fCBasket = GetNextBasket();
         return 1;
      }
   }
   return 0;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrack(GeantTrack &track, Bool_t priority)
{
// Add a track to the volume basket manager. If the track number reaches the
// threshold, the basket is added to the feeder queue and replaced by an empty 
// one. The feeder must be defined beforehand. Returns the number of dispatched
// baskets
#ifdef __STAT_DEBUG
   fScheduler->GetPendingStat().AddTrack(track);
#endif   
   if (priority) {
      fPBasket->AddTrack(track);
      if (fPBasket->GetNinput() >= fThreshold) {
#ifdef __STAT_DEBUG
         fScheduler->GetPendingStat().RemoveTracks(fPBasket->GetInputTracks());
         fScheduler->GetQueuedStat().AddTracks(fPBasket->GetInputTracks());
#endif   
         fFeeder->push(fPBasket, priority);
         fPBasket = GetNextBasket();
         return 1;
      }
   } else {
      fCBasket->AddTrack(track);
      if (fCBasket->GetNinput() >= fThreshold) {
#ifdef __STAT_DEBUG
         fScheduler->GetPendingStat().RemoveTracks(fCBasket->GetInputTracks());
         fScheduler->GetQueuedStat().AddTracks(fCBasket->GetInputTracks());
#endif   
         fFeeder->push(fCBasket, priority);
         fCBasket = GetNextBasket();
         return 1;
      }
   }
   return 0;
}      

//______________________________________________________________________________
Int_t GeantBasketMgr::CollectPrioritizedTracks(Int_t evmin, Int_t evmax)
{
// Move current basket tracks to priority one. 
// *** NONE *** This should be done for all basket managers only once when 
// starting prioritizing an event range.
   GeantTrack_v &tracks = fCBasket->GetInputTracks();
   Int_t ntracks = tracks.GetNtracks();
   for (Int_t itr=0; itr<ntracks; itr++) {
      if (tracks.fEventV[itr]>=evmin && tracks.fEventV[itr]<=evmax) {
         fPBasket->GetInputTracks().AddTracks(tracks, 0, ntracks-1);
         tracks.Clear();
         if (fPBasket->GetNinput() >= fThreshold) {
#ifdef __STAT_DEBUG
            fScheduler->GetPendingStat().RemoveTracks(fPBasket->GetInputTracks());
            fScheduler->GetQueuedStat().AddTracks(fPBasket->GetInputTracks());
#endif   
            fFeeder->push(fPBasket, kFALSE);
            fPBasket = GetNextBasket();
            return 1;
         }
         return 0;
      }
   }
   return 0;
}      

//______________________________________________________________________________
Int_t GeantBasketMgr::FlushPriorityBasket()
{
// Flush the baskets containing tracks. Returns the number of dispatched baskets.
   if (fPBasket->GetNinput()) {
#ifdef __STAT_DEBUG
            fScheduler->GetPendingStat().RemoveTracks(fPBasket->GetInputTracks());
            fScheduler->GetQueuedStat().AddTracks(fPBasket->GetInputTracks());
#endif   
      fFeeder->push(fPBasket, kTRUE);
      fPBasket = GetNextBasket();
      return 1;
   }
   return 0;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::GarbageCollect()
{
// Copy all priority tracks to the current basket and flush to queue
   GeantTrack_v &tracks = fPBasket->GetInputTracks();
   Int_t ntracks = tracks.GetNtracks();
   // Copy prioritized tracks to current basket
   if (ntracks) {
      fCBasket->GetInputTracks().AddTracks(tracks, 0, ntracks-1);
      tracks.Clear();
   }
   if (fCBasket->GetNinput()) {
#ifdef __STAT_DEBUG
      fScheduler->GetPendingStat().RemoveTracks(fCBasket->GetInputTracks());
      fScheduler->GetQueuedStat().AddTracks(fCBasket->GetInputTracks());
#endif   
      fFeeder->push(fCBasket, kFALSE);
      fCBasket = GetNextBasket();
      return 1;
   }
   return 0;
}   

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::GetNextBasket()
{
// Returns next empy basket if any available, else create a new basket.
   GeantBasket *next = fBaskets.try_pop();
   if (!next) {
      next = new GeantBasket(fThreshold+1, this);
      fMutex.Lock();
      // === critical section ===
      fNbaskets++;
      fNused++;
      // === end critical section ===
      fMutex.UnLock();
   } else {
      fMutex.Lock();
      // === critical section ===
      fNused++;
      // === end critical section ===
      fMutex.UnLock();
   }
   return next;
}

//______________________________________________________________________________
void GeantBasketMgr::RecycleBasket(GeantBasket *b)
{
// Recycle a basket.
   if (b->GetNinput() || b->GetNoutput()) {
      Printf("RecycleBasket: Error: ntracks!=0");
   }   
   b->Clear();
   if (b->GetInputTracks().Capacity() < fThreshold) {
      b->GetInputTracks().Resize(fThreshold);
      // Resize also the output array if needed
      if (b->GetOutputTracks().Capacity() < fThreshold)
        b->GetOutputTracks().Resize(fThreshold); 
   }     
   fBaskets.push(b);
   fMutex.Lock();
   // === critical section ===
   fNused--;
   // === end critical section ===
   fMutex.UnLock();
}   

//______________________________________________________________________________
void GeantBasketMgr::Print(Option_t *) const
{
// Print info about the basket content.
}

