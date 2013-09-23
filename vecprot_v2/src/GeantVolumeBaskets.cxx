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
                   fBaskets(),
                   fMutex()
{
// Constructor
}

//______________________________________________________________________________
GeantVolumeBaskets::~GeantVolumeBaskets()
{
// Clean up
   delete fBaskets;
}   

//______________________________________________________________________________
void GeantVolumeBaskets::GetNextBasket()
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

