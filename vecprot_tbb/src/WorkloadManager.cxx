#include "WorkloadManager.h"

#include "TGeoManager.h"

#include "GeantVolumeBasket.h"
#include "GeantMainScheduler.h"

WorkloadManager *WorkloadManager::fgInstance = 0;

//______________________________________________________________________________
WorkloadManager::WorkloadManager()
{
// Private constructor.
   fNbaskets = 0;
   fBasketArray = 0;
   AddEmptyCollections(32);
   fgInstance = this;
}

//______________________________________________________________________________
WorkloadManager::~WorkloadManager()
{
// Destructor.
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
   Int_t icrt = 0;
   while ((vol=(TGeoVolume*)next())) {
      basket = new GeantVolumeBasket(vol, icrt++);
      vol->SetField(basket);
      AddBasket(basket);
   }

   fMainSch = GeantMainScheduler::Instance(fNbaskets);
}

//______________________________________________________________________________
WorkloadManager *WorkloadManager::Instance()
{
// Return singleton instance.
   if (!fgInstance) fgInstance = new WorkloadManager();
   return fgInstance;
}

//______________________________________________________________________________
void WorkloadManager::AddEmptyBaskets(Int_t nb)
{
// Add empty baskets in the queue.
   for (Int_t i=0; i<nb; i++) tbb_feeder_empty_queue.push(new GeantBasket(10));
   Printf("Added %d empty baskets to the queue", nb);
}

//______________________________________________________________________________
void WorkloadManager::AddEmptyCollections(Int_t nc)
{
// Add empty baskets in the queue.
   for (Int_t i=0; i<nc; i++) tbb_collector_empty_queue.push(new GeantTrackCollection(100));
   Printf("Added %d empty collections to the queue", nc);
}


