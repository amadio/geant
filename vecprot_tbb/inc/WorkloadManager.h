#ifndef GEANT_WORKLOADMANAGER
#define GEANT_WORKLOADMANAGER

#include "Rtypes.h"

#include "tbb/concurrent_queue.h"

#include "GeantBasket.h"

class GeantVolumeBasket;
class GeantMainScheduler;

// Main work manager class. This creates and manages all the worker threads,
// has pointers to the synchronization objects, but also to the currently
// transported baskets.

//______________________________________________________________________________
class WorkloadManager
{
public:
   Int_t              fNthreads;           // number of managed threads
   Int_t              fNbaskets;           // total number of baskets
   Int_t              fNbasketgen;         // number of baskets to transport in the current generation
   Int_t              fNidle;              // number of idle workers
   Int_t              fNqueued;            // Number of chunks queued
   Int_t             *fBtogo;              // array of baskets to be processed in the next generation
   Bool_t             fStarted;            // Start flag

   tbb::concurrent_bounded_queue<GeantBasket*> tbb_feeder_queue;
   tbb::concurrent_bounded_queue<GeantBasket*> tbb_feeder_priority_queue;
   tbb::concurrent_bounded_queue<GeantBasket*> tbb_feeder_empty_queue;
   tbb::concurrent_bounded_queue<GeantTrackCollection*> tbb_collector_queue;
   tbb::concurrent_bounded_queue<GeantTrackCollection*> tbb_collector_empty_queue;

   static WorkloadManager *fgInstance;     // Singleton instance
   GeantVolumeBasket **fBasketArray;       //![number of volumes] Array of baskets
   Bool_t             fFlushed;            // Buffer flushed
   Bool_t             fFilling;            // Worker queue is filling

public:
   WorkloadManager(Int_t nthreads);
   virtual ~WorkloadManager();
   void                AddBasket(GeantVolumeBasket *basket) {fBasketArray[fNbaskets++]=basket;}
   void                AddEmptyBaskets(Int_t nb);
   void                CreateBaskets(Int_t nvolumes);

   Int_t               GetNthreads() const {return fNthreads;}
   Int_t               GetNbaskets() const {return fNbaskets;}
   GeantVolumeBasket **GetBasketArray() const {return fBasketArray;}
   static WorkloadManager *
                       Instance(Int_t nthreads=0);

   GeantMainScheduler* fMainSch;

   Bool_t              IsFlushed() const {return fFlushed;}
   Bool_t              IsFilling() const {return fFilling;}
   void                SetFlushed(Bool_t flag) {fFlushed = flag;}
};
#endif // GEANT_WORKLOADMANAGER

