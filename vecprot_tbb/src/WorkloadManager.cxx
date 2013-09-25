#include "WorkloadManager.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TBits.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TThread.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TGeoBranchArray.h"

#include "GeantPropagator.h"
#include "GeantVolumeBasket.h"
#include "GeantOutput.h"
#include "GeantThreadData.h"
#include "PhysicsProcess.h"
#include "Event_Hit.h"
#include "GeantMainScheduler.h"

WorkloadManager *WorkloadManager::fgInstance = 0;

//______________________________________________________________________________
WorkloadManager::WorkloadManager(Int_t nthreads)
{
// Private constructor.
   fNthreads = nthreads;
   fNbaskets = 0;
   fNbasketgen = 0;
   fNidle = nthreads;
   fNqueued = 0;
   fBtogo = 0;
   fStarted = kFALSE;

   for (Int_t i=0; i<nthreads; i++) tbb_collector_empty_queue.push(new GeantTrackCollection(100));

   fStarted = kFALSE;
   fgInstance = this;
   fBasketArray = 0;
   fFlushed = kFALSE;
   fFilling = kFALSE;
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
void WorkloadManager::AddEmptyBaskets(Int_t nb)
{
// Add empty baskets in the queue.
   for (Int_t i=0; i<nb; i++) tbb_feeder_empty_queue.push(new GeantBasket(10));
   Printf("Added %d empty baskets to the queue", nb);
}

