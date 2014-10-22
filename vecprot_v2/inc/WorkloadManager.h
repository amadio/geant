#ifndef GEANT_WORKLOADMANAGER
#define GEANT_WORKLOADMANAGER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "priority_queue.h"
#include "rr_pool.h"
#include "GeantObjectPool.h"
#include "GeantTrack.h"
 
class GeantBasketMgr;
class GeantBasket;
class GeantScheduler;
class TaskBroker;

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
   Bool_t             fStopped;            // Stop flag
   Bool_t             fWorkDone;           // Flag that a worker has published some result
   priority_queue<GeantBasket*> 
                     *fFeederQ;            // queue of transportable baskets
   priority_queue<GeantBasket*>  
                     *fTransportedQ;       // queue of transported baskets
   priority_queue<GeantBasket*> 
                     *fDoneQ;              // Thread "all work done" queue
//   GeantObjectPool<VolumePath_t>
//                     *fNavStates;          // Pool of navigation states                  
   rr_pool<VolumePath_t>
                     *fNavStates;          // Pool of navigation states                  
   static WorkloadManager *fgInstance;     // Singleton instance
   TList             *fListThreads;        // List of threads
   Bool_t             fFlushed;            // Buffer flushed
   Bool_t             fFilling;            // Worker queue is filling
   GeantScheduler    *fScheduler;          // Main basket scheduler

   TaskBroker        *fBroker;             // Pointer to the coprocessor broker, this could be made a collection.
   Int_t             *fWaiting;            //![fNthreads+1] Threads in waiting flag
#if __cplusplus >= 201103L
   std::mutex        *fMutexSch;           // mutex for the scheduler
   std::condition_variable
                     *fCondSch;            // Wait condition for scheduler
#endif
   WorkloadManager(Int_t nthreads);
public:
   virtual ~WorkloadManager();
   void                CreateBaskets();
   priority_queue<GeantBasket*> *FeederQueue() const {return fFeederQ;}
   priority_queue<GeantBasket*> *TransportedQueue() const {return fTransportedQ;}
   priority_queue<GeantBasket*> *DoneQueue() const {return fDoneQ;}
//   GeantObjectPool<VolumePath_t>  
   rr_pool<VolumePath_t>  
                      *NavStates() const   {return fNavStates;}
   Int_t               GetNthreads() const {return fNthreads;}
   Int_t               GetNbaskets() const {return fNbaskets;}
   Int_t              *GetWaiting() const  {return fWaiting;}
   GeantScheduler     *GetScheduler() const {return fScheduler;}

#if __cplusplus >= 201103L
   std::mutex         *GetMutexSch() const {return fMutexSch;}
   std::condition_variable *GetCondSch() const {return fCondSch;}
#endif
   static WorkloadManager *
                       Instance(Int_t nthreads=0);                    
   Bool_t              IsFlushed() const {return fFlushed;}
   Bool_t              IsFilling() const {return fFilling;}
   Bool_t              IsStopped() const {return fStopped;}
   Bool_t              IsWorkDone() const {return fWorkDone;}
   void                SetWorkDone(Bool_t flag) {fWorkDone = flag;}
   void                Stop()            {fStopped = kTRUE;}
   void                SetFlushed(Bool_t flag) {fFlushed = flag;}
   Int_t               GetBasketGeneration() const {return fBasketGeneration;}
   void                Print(Option_t *option="") const;
   void                SetTaskBroker(TaskBroker *broker);
   Int_t               GetNminThreshold() const {return fNminThreshold;}
   void                SetNminThreshold(Int_t thr) {fNminThreshold = thr;}
   void                StartThreads();
   void                JoinThreads();
   static void        *MainScheduler(void *arg);
   static void        *MonitoringThread(void *arg);
   static void        *TransportTracks(void *arg);
   static void        *TransportTracksCoprocessor(void *arg);
   void                WaitWorkers();
private:
   WorkloadManager(const WorkloadManager &);//no imp.	
   WorkloadManager& operator=(const WorkloadManager &);//no imp.
      
   
   ClassDef(WorkloadManager,0)  // The work manager class.
};   
#endif
