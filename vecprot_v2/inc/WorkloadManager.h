#ifndef GEANT_WORKLOADMANAGER
#define GEANT_WORKLOADMANAGER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "sync_objects.h"
 
class GeantBasketMgr;
class GeantBasket;
class GeantScheduler;

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
   dcqueue<GeantBasket> 
                     *fFeederQ;            // queue of transportable baskets
   dcqueue<GeantBasket>  
                     *fTransportedQ;       // queue of transported baskets
   dcqueue<GeantBasket> 
                     *fDoneQ;              // Thread "all work done" queue
   static WorkloadManager *fgInstance;     // Singleton instance
   TList             *fListThreads;        // List of threads
   Bool_t             fFlushed;            // Buffer flushed
   Bool_t             fFilling;            // Worker queue is filling
   GeantScheduler    *fScheduler;          // Main basket scheduler

   WorkloadManager(Int_t nthreads);
public:
   virtual ~WorkloadManager();
   void                CreateBaskets();
   dcqueue<GeantBasket> *FeederQueue() const {return fFeederQ;}
   dcqueue<GeantBasket> *TransportedQueue() const {return fTransportedQ;}
   dcqueue<GeantBasket> *DoneQueue() const {return fDoneQ;}
   
   Int_t               GetNthreads() const {return fNthreads;}
   Int_t               GetNbaskets() const {return fNbaskets;}
   GeantScheduler     *GetScheduler() const {return fScheduler;}
   static WorkloadManager *
                       Instance(Int_t nthreads=0);                    
   Bool_t              IsFlushed() const {return fFlushed;}
   Bool_t              IsFilling() const {return fFilling;}
   void                SetFlushed(Bool_t flag) {fFlushed = flag;}
   Int_t               GetBasketGeneration() const {return fBasketGeneration;}
   void                Print(Option_t *option="") const;
   Int_t               GetNminThreshold() const {return fNminThreshold;}
   void                SetNminThreshold(Int_t thr) {fNminThreshold = thr;}
   void                StartThreads();
   void                JoinThreads();
   static void        *MainScheduler(void *arg);
   static void        *TransportTracks(void *arg);
   static void        *TransportTracksCoprocessor(void *arg);
   void                WaitWorkers();
   
   ClassDef(WorkloadManager,0)  // The work manager class.
};   
#endif
