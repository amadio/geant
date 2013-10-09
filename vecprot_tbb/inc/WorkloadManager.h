#ifndef GEANT_WORKLOADMANAGER
#define GEANT_WORKLOADMANAGER

#include "Rtypes.h"

#include "tbb/concurrent_queue.h"

#include "GeantBasket.h"

class GeantVolumeBasket;
class GeantMainScheduler;

//______________________________________________________________________________
class WorkloadManager
{
public:
   static WorkloadManager *fgInstance;

   Int_t fNbaskets;
   GeantVolumeBasket **fBasketArray;

   tbb::concurrent_bounded_queue<GeantBasket*> tbb_feeder_queue;
   tbb::concurrent_bounded_queue<GeantBasket*> tbb_feeder_priority_queue;
   tbb::concurrent_bounded_queue<GeantBasket*> tbb_feeder_empty_queue;
   tbb::concurrent_bounded_queue<GeantTrackCollection*> tbb_collector_queue;
   tbb::concurrent_bounded_queue<GeantTrackCollection*> tbb_collector_empty_queue;

public:
   WorkloadManager();
   virtual ~WorkloadManager();
   void AddBasket(GeantVolumeBasket *basket) {fBasketArray[fNbaskets++]=basket;}
   void AddEmptyBaskets(Int_t nb);
   void AddEmptyCollections(Int_t nc);
   void CreateBaskets(Int_t nvolumes);

   Int_t GetNbaskets() const {return fNbaskets;}
   GeantVolumeBasket **GetBasketArray() const {return fBasketArray;}
   static WorkloadManager* Instance();

   GeantMainScheduler* fMainSch;

};

#endif // GEANT_WORKLOADMANAGER

