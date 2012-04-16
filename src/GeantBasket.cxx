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
void GeantBasket::Print(Option_t *) const
{
// Print basket content.
   for (Int_t i=0; i<fNtracks; i++) {
      gPropagator->fTracks[fIndex[i]]->Print();
   }
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
   
ClassImp(GeantBasketScheduler)

//______________________________________________________________________________
GeantBasketScheduler::GeantBasketScheduler()
        :TObject(),
         fVolume(0),
         fNbaskets(0),
         fMaxBaskets(0),
         fVIB(0),
         fCrt(0),
         fNbusy(0),
         fBaskets(0),
         fCollector(0)
{
// ctor
}

//______________________________________________________________________________
GeantBasketScheduler::GeantBasketScheduler(TGeoVolume *vol)
        :TObject(),
         fVolume(vol),
         fNbaskets(0),
         fMaxBaskets(10),
         fVIB(0),
         fCrt(0),
         fNbusy(0),
         fBaskets(0),
         fCollector(0)
{
// ctor
   fBaskets = new GeantBasket*[fMaxBaskets];
   for (Int_t i=0; i<fMaxBaskets-1; i++) AddNewBasket();
   if (gPropagator->fWMgr->EmptyQueue()->empty())
      fCollector = new GeantBasket(1000);
   else { 
      fCollector = gPropagator->fWMgr->EmptyQueue()->wait_and_pop();
      fCollector->Clear();
   }   
}

//______________________________________________________________________________
GeantBasketScheduler::~GeantBasketScheduler()
{
// Destructor
   if (fBaskets) {
      for (Int_t i=0; i<fNbaskets; i++) delete fBaskets[i];
   }
   delete [] fBaskets;
   delete fCollector;
}   

//______________________________________________________________________________
void GeantBasketScheduler::AddNewBasket()
{
// Add new basket in the round robin. Pick one from the queue of empty baskets
// if available, otherwise create new one.
   if (fNbaskets >= fMaxBaskets-1) {
      GeantBasket **newarray = new GeantBasket*[2*fMaxBaskets];
      memcpy(newarray, fBaskets, fNbaskets*sizeof(GeantBasket*));
      delete [] fBaskets;
      fBaskets = newarray;
      fMaxBaskets *= 2;
   }
   if (gPropagator->fWMgr->EmptyQueue()->empty())
      fBaskets[fNbaskets++] = new GeantBasket(10);
   else { 
      fBaskets[fNbaskets++] = gPropagator->fWMgr->EmptyQueue()->wait_and_pop();
      fBaskets[fNbaskets-1]->Clear();
   }   
}

//______________________________________________________________________________
void GeantBasketScheduler::AddTrack(Int_t itrack)
{
// Preferential round robin. The prefered basket is pointed by fVIB. The first
// thread comung will find fCrt = fVIB. It will atomically increment the round 
// robin fCrt for the next thread, while filling the track to fVIB itself. The
// next concurrent thread will write in the next slot and so on. After writing
// the track to the basket, the thread decremnts fNbusy and changes fCrt = fVIB
// if it is the last worker on this basket.

   // Read the current value of fCrt
   Int_t crt = 0;
   fBasketLock.Lock();
/**/ fNbusy++;
/**/ crt = fCrt;
/**/ // The line below should not happen too often...only for many threads and in 
/**/ // the training phase => TO BE CHECKED
/**/ if (fNbusy>fNbaskets) AddNewBasket();
/**/ // Round robin the basket to be filled
/**/ fCrt = (fCrt+1)%fNbaskets;
   if (fCrt==crt) printf("fCrt=%d fNbusy=%d\n",fCrt, fNbusy);
   fBasketLock.UnLock(); // no other thread should work with icrt after this point
   // Add track to slot pointed by local variable crt
   fBaskets[crt]->AddTrack(itrack);
   // Check if basket became transportable
   if (fBaskets[crt]->GetNtracks() >= gPropagator->fNperBasket) {
      InjectBasket(crt); // only one writing there
      if (fVIB==crt) {
         fBasketLock.Lock();
         // Round robin the priority basket
/**/     fVIB = (fVIB+1)%fNbaskets;
         fBasketLock.UnLock();
      }   
   }
   fBasketLock.Lock();
/**/ fNbusy--;
   // If last in line, make the VIB current
/**/ if (!fNbusy) fCrt = fVIB;
   fBasketLock.UnLock();
}   

//______________________________________________________________________________
GeantBasket *GeantBasketScheduler::GarbageCollect(GeantBasket *basket, Bool_t feed)
{
// Garbage collect baskets managed by this scheduler. This method is called
// ONLY by the garbage collector thread. An empty basket is provided as input.
// Garbage collect only if the scheduler is not busy and has tracks. Returns
// recycled basket;
   if (fNbusy) return basket;
   if (feed && !GetNtracks()) {
      if (!fCollector->GetNtracks()) return basket;
      gPropagator->fWMgr->FeederQueue()->push(fCollector);
      fCollector = basket;
      return 0;
   }   
      
   GeantBasket *btmp = basket;
   Int_t ntracks = 0;
   // Loop baskets starting from (fVIB-1)%fNbaskets backwards
   for (Int_t i=0; i<fNbaskets; i++) {
      Int_t ib = (fVIB+fNbaskets-i-1)%fNbaskets;
      if (!fBaskets[ib]->GetNtracks()) continue;
      fBasketLock.Lock();
/**/  if (fNbusy) {
         fBasketLock.UnLock();
         return btmp;
      }
      // Steal the basket and replace with empty one
/**/  btmp = fBaskets[ib];
/**/  ntracks = btmp->GetNtracks();
/**/  fBaskets[ib] = basket;
      fBasketLock.UnLock();
      // btmp points to basket to be flushed and recycled
      fCollector->AddTracks(btmp->GetTracks(), ntracks);
      btmp->Clear();
   }
   // Inject collector if not empty and replace with the recycled one.
   if (!feed || !fCollector->GetNtracks()) return btmp;
   gPropagator->fWMgr->FeederQueue()->push(fCollector);
   fCollector = btmp;
   return 0;
}   

//______________________________________________________________________________
Int_t GeantBasketScheduler::GetNtracks() const
{
// Get a snapshot of the dynamically changing number of tracks. Note that
// the number of tracks can change DURING the call. This is OK since this
// number is only used to trigger basket injection.
   Int_t ntracks = 0;
   Int_t nb = fNbaskets; // the number may only increase
   for (Int_t ib=0; ib<nb; ib++) ntracks += fBaskets[ib]->GetNtracks();
   return ntracks;
}

//______________________________________________________________________________
Int_t GeantBasketScheduler::GetNtotal() const
{
   return GetNtracks()+fCollector->GetNtracks();
}
   
//______________________________________________________________________________
void GeantBasketScheduler::InjectBasket(Int_t islot)
{
// Inject basket at islot in the processing queue and replace with an empty one.
// TO CHECK: can the push operation become blocking ?
   GeantBasket *empty = 0;
   if (gPropagator->fWMgr->EmptyQueue()->empty())
      empty = new GeantBasket(10);
   else { 
      empty = gPropagator->fWMgr->EmptyQueue()->wait_and_pop();
      empty->Clear();
   }   
   GeantBasket *full = fBaskets[islot];
   fBasketLock.Lock();
   fBaskets[islot] = empty;
   fBasketLock.UnLock();
   gPropagator->fWMgr->FeederQueue()->push(full);
}

