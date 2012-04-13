#ifndef GEANT_WORKLOADMANAGER
#define GEANT_WORKLOADMANAGER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "sync_objects.h"
 
class GeantVolumeBasket;

// Main work manager class. This creates and manages all the worker threads,
// has pointers to the synchronization objects, but also to the currently
// transported baskets.

//______________________________________________________________________________
class WorkloadManager : public TObject {
protected:
   Int_t              fNthreads;           // number of managed threads
   Int_t              fNbaskets;           // total number of baskets
   Int_t              fBasketGeneration;   // basket generation
   Int_t              fNbasketgen;         // number of baskets to transport in the current generation
   Int_t              fNidle;              // number of idle workers
   Int_t              fNminThreshold;      // Minimum number of tracks in a basket to trigger transport
   Int_t              fNqueued;            // Number of chunks queued
   Int_t             *fBtogo;              // array of baskets to be processed in the next generation
   Bool_t             fStarted;            // Start flag
   concurrent_queue  *feeder_queue;        // queue of filled baskets
   concurrent_queue  *answer_queue;        // queue of filled baskets
   concurrent_queue  *empty_queue;         // queue of empty baskets
   static WorkloadManager *fgInstance;     // Singleton instance
   GeantVolumeBasket **fCurrentBasket;     // Current basket transported per thread
   TList             *fListThreads;        // List of threads
   GeantVolumeBasket **fBasketArray;       //![number of volumes] Array of baskets
   Bool_t             fFlushed;            // Buffer flushed
   Bool_t             fFilling;            // Worker queue is filling

   WorkloadManager(Int_t nthreads);
public:
   virtual ~WorkloadManager();
   void                AddBasket(GeantVolumeBasket *basket) {fBasketArray[fNbaskets++]=basket;}
   void                AddEmptyBaskets(Int_t nb);
   void                InterruptBasket(GeantVolumeBasket *basket, Int_t *trackin, Int_t ntracks, Int_t tid);
   void                CreateBaskets(Int_t nvolumes);
   concurrent_queue   *AnswerQueue() const {return answer_queue;}
   concurrent_queue   *EmptyQueue() const  {return empty_queue;}
   concurrent_queue   *FeederQueue() const {return feeder_queue;}
   
   Int_t               GetNthreads() const {return fNthreads;}
   Int_t               GetNbaskets() const {return fNbaskets;}
   GeantVolumeBasket **GetBasketArray() const {return fBasketArray;}
   static WorkloadManager *
                       Instance(Int_t nthreads=0);
   Bool_t              IsFlushed() const {return fFlushed;}
   Bool_t              IsFilling() const {return fFilling;}
   void                SetFlushed(Bool_t flag) {fFlushed = flag;}
   Int_t               GetBasketGeneration() const {return fBasketGeneration;}
   void                SetCurrentBasket(Int_t tid, GeantVolumeBasket *basket) {fCurrentBasket[tid]=basket;}
   GeantVolumeBasket  *GetCurrentBasket(Int_t tid) const {return fCurrentBasket[tid];}
   void                Print();
   void                SetNminThreshold(Int_t thr) {fNminThreshold = thr;}
   void                StartThreads();
   void                JoinThreads();
   static void        *GarbageCollect(void *arg);
   static void        *TransportTracks(void *arg);
   void                WaitWorkers();
   
   ClassDef(WorkloadManager,0)  // The work manager class.
};   
#endif
