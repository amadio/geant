#ifndef GEANT_WORKLOADMANAGER
#define GEANT_WORKLOADMANAGER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "sync_objects.h"
 
class GeantVolumeBasket;
class GeantParticleBuffer;

// Main work manager class. This creates and manages all the worker threads,
// has pointers to the synchronization objects, but also to the currently
// transported baskets.

//______________________________________________________________________________
class WorkloadManager : public TObject {
protected:
   Int_t              fNthreads;              // number of managed threads
   Int_t              fNbaskets;              // total number of baskets
   Int_t              fBasketGeneration;      // basket generation
   Int_t              fNbasketgen;            // number of baskets to transport in the current generation
   Int_t              fNidle;                 // number of idle workers
   Int_t              fNminThreshold;         // Minimum number of tracks in a basket to trigger transport
   Int_t              fNqueued;               // Number of chunks queued
   Int_t             *fBindex;                // array of basket indices in the current generation
   Bool_t             fStarted;               // Start flag
   concurrent_queue  *feeder_queue;           // workload queue
   concurrent_queue  *answer_queue;           // answer queue
   static WorkloadManager *fgInstance;        // Singleton instance
   GeantVolumeBasket **fCurrentBasket;        // Current basket transported per thread
   TList             *fListThreads;          // List of threads
   GeantVolumeBasket **fBasketArray; //![number of volumes] Array of baskets
   GeantParticleBuffer *fBuffer;     // Buffer of pending particles

   WorkloadManager(Int_t nthreads);
public:
   virtual ~WorkloadManager();
   void                AddBasket(GeantVolumeBasket *basket) {fBasketArray[fNbaskets++]=basket;}
   void                AddPendingTrack(Int_t itrack, GeantVolumeBasket *basket, Int_t tid);
   void                CreateBaskets(Int_t nvolumes);
   concurrent_queue   *FeederQueue() const {return feeder_queue;}
   concurrent_queue   *AnswerQueue() const {return answer_queue;}
   
   static WorkloadManager *
                       Instance(Int_t nthreads=0);
   void                ClearBaskets();
   void                QueueBaskets();
   Int_t               GetBasketGeneration() const {return fBasketGeneration;}
   GeantParticleBuffer *GetBuffer() const {return fBuffer;}
   void                SetCurrentBasket(Int_t tid, GeantVolumeBasket *basket) {fCurrentBasket[tid]=basket;}
   GeantVolumeBasket  *GetCurrentBasket(Int_t tid) const {return fCurrentBasket[tid];}
   void                Print();
   void                SortBaskets();
   void                SelectBaskets();
   void                StartThreads();
   void                JoinThreads();
   static void        *TransportTracks(void *arg);
   void                WaitWorkers();
   
   ClassDef(WorkloadManager,0)  // The work manager class.
};   
#endif
