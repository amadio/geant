#include "GeantVolumeBaskets.h"
#include "GeantTrack.h"

#include "TThread.h"

//______________________________________________________________________________
// Basket manager for a given volume. Holds a list of free baskets stored in a
// concurrent queue
//______________________________________________________________________________


ClassImp(GeantVolumeBaskets)

const Double_t gTolerance = TGeoShape::Tolerance();

//______________________________________________________________________________
GeantVolumeBaskets::GeantVolumeBasket(TGeoVolume *vol, Int_t number)
                  :TObject(),
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
}

//______________________________________________________________________________
GeantVolumeBaskets::~GeantVolumeBaskets()
{
// Clean up
   delete fBaskets;
   delete fCBasket;
   delete fPBasket;
}   

//______________________________________________________________________________
Int_t GeantVolumeBaskets::AddTrack(const GeantTrack &track, Bool_t priority)
{
// Add a track to the volume basket manager. If the track number reaches the
// threshold, the basket is added to the feeder queue and replaced by an empty 
// one. The feeder must be defined beforehand. Returns the number of dispatched
// baskets
   GeantBasket *basket = 0;
   if (priority) {
      if (!fPBasket) fPBasket = GetNextBasket();
      fPBasket->AddTrack(track);
      if (fPBasket->GetNinput() >= fThreshold) {
         fFeeder->push(fPBasket, priority);
         fPBasket = GetNextBasket();
         return 1;
      }
   } else {
      if (!fCBasket) fCBasket = GetNextBasket();
      fCBasket->AddTrack(track);
      if (fCBasket->GetNinput() >= fThreshold) {
         fFeeder->push(fCBasket, priority);
         fCBasket = GetNextBasket();
         return 1;
      }
   }
   return 0;
}      

//______________________________________________________________________________
Int_t GeantVolumeBaskets::Flush()
{
// Flush the baskets containing tracks. Returns the number of dispatched baskets.
   Int_t ndispatched = 0;
   if (fPBasket && fPBasket->GetNinput()) {
      fFeeder->push(fPBasket, kTRUE);
      fPBasket = GetNextBasket();
      ndispatched++;
   }
   if (fCBasket && fCBasket->GetNinput()) {
      fFeeder->push(fCBasket, kFALSE);
      fCBasket = GetNextBasket();
      ndispatched++;
   }
   return ndispatched;
}

//______________________________________________________________________________
GeantBasket *GeantVolumeBaskets::GetNextBasket()
{
// Returns next empy basket if any available, else create a new basket.
   GeantBasket *next = fBaskets.try_pop();
   if (!next) {
      next = new GeantBasket(fThreshold+1, fNumber, 
                 fVolume->GetMaterial()->GetNumber()); // maybe bigger, don't know...
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
void GeantVolumeBaskets::RecycleBasket(GeantBasket *b)
{
// Recycle a basket.
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
void GeantVolumeBaskets::Print(Option_t *) const
{
// Print info about the basket content.
}

