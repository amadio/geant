//===--- WorkloadManager.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file WorkloadManager.h
 * @brief Definition of workload manager in Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_WORKLOADMANAGER
#define GEANT_WORKLOADMANAGER

#include <thread>
#include <list>
#include "priority_queue.h"
#include "dcqueue.h"
#include "condition_locker.h"

#ifndef USE_VECGEOM_NAVIGATOR
#include "TGeoManager.h"
#endif

#include "GeantTrackVec.h"
#include "GeantPropagator.h"

#include "GeantTaskData.h"

#ifdef USE_ROOT
#include "TThreadMergingServer.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantBasketMgr;
class GeantBasket;
class GeantScheduler;
class TaskBroker;
class GeantVTaskMgr;
class EventSet;

/**
 * @brief WorkloadManager class
 * @details Main work manager class. This creates and manages all the worker threads,
 * has pointers to the synchronization objects, but also to the currently
 * transported baskets.
 */
class WorkloadManager {
protected:
  GeantPropagator* fPropagator;                           /** propagator,... */
  int fNthreads;                                     /** Number of managed threads */
  int fNbaskets;                                     /** Total number of baskets */
  int fBasketGeneration;                             /** Basket generation */
  int fNbasketgen;                                   /** Number of baskets to transport in the current generation */
  int fNidle;                                        /** Number of idle workers */
  int fNqueued;                                      /** Number of chunks queued */
  int *fBtogo;                                       /** Array of baskets to be processed in the next generation */
  int fSchId;                                        /** Thread id for the scheduler */
  bool fStarted;                                       /** Start flag */
  bool fStopped;                                       /** Stop flag */
  Geant::priority_queue<GeantBasket *> *fFeederQ;      /** Queue of transportable baskets */
  Geant::priority_queue<GeantBasket *> *fTransportedQ; /** Queue of transported baskets */
  Geant::priority_queue<GeantBasket *> *fDoneQ;        /** Thread "all work done" queue */
  std::vector<std::thread> fListThreads;               /** Vector of threads */
  bool fFlushed;                                       /** Buffer flushed */
  bool fFilling;                                       /** Worker queue is filling */
  GeantScheduler *fScheduler;                          /** Main basket scheduler */

  TaskBroker *fBroker;         /** Pointer to the coprocessor broker, this could be made a collection. */
  condition_locker fSchLocker; /** Scheduler locker */
  condition_locker fGbcLocker; /** Garbage collector locker */
  std::atomic_flag fShareLock = ATOMIC_FLAG_INIT; /** Atomic flag to protect basket sharing */
  int fLastEvent;            /** Last transported event */

  #ifdef USE_ROOT
  dcqueue<TBufferFile*>* fOutputIO;            /** Queue of buffers to be merged for IO **/ 
  Geant::TThreadMergingServer* fMergingServer;
  #endif

  /**
   * @brief WorkloadManager parameterized constructor
   *
   * @param  nthreads Number of threads foe workload manager
   */
  WorkloadManager(int nthreads, GeantPropagator* prop);

public:
  /** @brief WorkloadManager destructor */
  virtual ~WorkloadManager();

  /** @brief Create basket function */
  void CreateBaskets(GeantPropagator* prop);

  enum class FeederResult : char { kNone, kError, kWork, kStop };

  /** @brief Call Feeder (if needed) and check exit condition. */
  FeederResult CheckFeederAndExit();

  /** @brief Function for feeder queue of transportable baskets */
  Geant::priority_queue<GeantBasket *> *FeederQueue() const { return fFeederQ; }

  /**
   * @brief Function that provides trasported queue for baskets
   */
  Geant::priority_queue<GeantBasket *> *TransportedQueue() const { return fTransportedQ; }

  /**
   * @brief Function that provides thread's "all work done" queue
   */
  Geant::priority_queue<GeantBasket *> *DoneQueue() const { return fDoneQ; }
  //   GeantObjectPool<VolumePath_t>
  //   rr_pool<VolumePath_t>
  //                      *NavStates() const   {return fNavStates;}

  /**
   * @brief Function that provides IO queue
   */
  #ifdef USE_ROOT
  dcqueue<TBufferFile*>*  IOQueue() const { return fOutputIO; }
  Geant::TThreadMergingServer* MergingServer() const { return fMergingServer; }
  #endif
  /** @brief Function that returns number of managed threads */
  GEANT_FORCE_INLINE
  int GetNthreads() const { return fNthreads; }

  /** @brief Function that returns total number of baskets */
  GEANT_FORCE_INLINE
  int GetNbaskets() const { return fNbaskets; }

  /** @brief Function that returns number of pending baskets */
  GEANT_FORCE_INLINE
  int GetNpending() const { return fFeederQ->size_async(); }

  ///** @brief Function returning the number of monitored features */
  //int GetMonFeatures() const;

  /** @brief Function that returns main basket scheduler */
  GEANT_FORCE_INLINE
  GeantScheduler *GetScheduler() const { return fScheduler; }

  /** @brief Get the scheduler thread id */
  GEANT_FORCE_INLINE
  int GetSchId() const { return fSchId; }

  /** @brief Set scheduler thread id */
  GEANT_FORCE_INLINE
  void SetSchId(int id) { fSchId = id; }

  /** @brief Function that returns scheduler locker */
  GEANT_FORCE_INLINE
  condition_locker &GetSchLocker() { return fSchLocker; }

  /** @brief Function that returns garbage collector locker */
  GEANT_FORCE_INLINE
  condition_locker &GetGbcLocker() { return fGbcLocker; }

  /**
   * @brief Function that create workload manager instance
   *
   * @param nthreads Number of threads (by default 0)
   */
  static WorkloadManager *NewInstance(GeantPropagator *prop= nullptr, int nthreads = 0);

  /** @brief Function that check if buffer is flushed */
  GEANT_FORCE_INLINE
  bool IsFlushed() const { return fFlushed; }

  /** @brief Function that check if worker queue is filling */
  GEANT_FORCE_INLINE
  bool IsFilling() const { return fFilling; }

  /** @brief Function that check stop flag */
  GEANT_FORCE_INLINE
  bool IsStopped() const { return fStopped; }

  /** @brief Getter for last transported event */
  GEANT_FORCE_INLINE
  int LastEvent() const { return fLastEvent; }

  /** @brief Setter for last transported event */
  GEANT_FORCE_INLINE
  void SetLastEvent(int n) { fLastEvent = n; }

  /** @brief Function that provide stop process by setting Stop flag = True */
  GEANT_FORCE_INLINE
  void Stop() { fStopped = true; }

  /** @brief Stop all transport threads */
  void StopTransportThreads();

  /** @brief Setter for buffer flushing */
  GEANT_FORCE_INLINE
  void SetFlushed(bool flag) { fFlushed = flag; }

  /** @brief Function that returns basket generation */
  GEANT_FORCE_INLINE
  int GetBasketGeneration() const { return fBasketGeneration; }

  /** @brief Print function */
  void Print(const char *option = "") const;

  /** @brief  Setter for task broker */
  void SetTaskBroker(TaskBroker *broker);

  /** @brief  Setter for task broker */
  GEANT_FORCE_INLINE
  TaskBroker *GetTaskBroker() { return fBroker; }

  /** @brief Try to acquire the sharing lock */
  GEANT_FORCE_INLINE
  bool TryShareLock() { return (fShareLock.test_and_set(std::memory_order_acquire)); }

  /** @brief Release the share lock */
  GEANT_FORCE_INLINE
  void ReleaseShareLock() { fShareLock.clear(std::memory_order_release); }

  /** @brief Check if transport is feeding with new tracks. */
  GEANT_FORCE_INLINE
  bool IsSharing() {
    if (TryShareLock()) return true;
    ReleaseShareLock();
    return false;
  }

#if USE_VECGEOM_NAVIGATOR
  /**
   * @brief Tell the task broker(s) to load the geometry.
   *
   * @param Volume to load
   */
  bool LoadGeometry(vecgeom::VPlacedVolume const *const volume = nullptr);
#endif

//   void SetMaxThreads(int nthreads) {
//   fMaxThreads = nthreads;
// #ifndef USE_VECGEOM_NAVIGATOR
//     gGeoManager->SetMaxThreads(nthreads);
// #endif
//   }

  int ThreadId();

  // /** @brief Getter for the global transport threshold */
  // int GetNminThreshold() const { return fNminThreshold; }

  // /** @brief Setter for the global transport threshold */
  // void SetNminThreshold(int thr) { fNminThreshold = thr; }

  /** @brief Function that initializes the threads/tasks used by the system */
  bool StartTasks(GeantVTaskMgr *taskmgr);

  /** @brief Joins all threads at the end of processing */
  void JoinThreads();

  /** @brief Thread function for the main scheduler */
  static void *MainScheduler(void *arg);

  /**
   * @brief Function that provides garbage collector thread
   *
   * @param arg Arguments to be passed in the function
   */
  static void *GarbageCollectorThread(GeantPropagator *prop);

  /**
   * @brief Function that provides monitoring thread
   *
   * @param arg Arguments to be passed in the function
   */
  static void *MonitoringThread(GeantPropagator *prop);

  /**
   * @brief Function that provides output thread
   *
   * @param arg Arguments to be passed in the function
   */  
  static void *OutputThread(GeantPropagator *prop);

#ifdef USE_ROOT
  /** @brief Function that starts ROOT application */  
  static void *StartROOTApplication();
#endif

  /**
   * @brief Function that provides transporting tracks
   *
   * @param arg Arguments to be passed in the function
   */
  static void *TransportTracks(GeantPropagator *prop);

  /**
   * @brief Function that provides transporting tracks
   *
   * @param arg Arguments to be passed in the function
   */
  static void TransportTracksV3(GeantPropagator *prop);

  /**
   * @brief Function that provides transporting tracks
   *
   * @param td Task data
   * @param workload Event set to be transported
   * @return True if workload completed. If false, the work will be completed
   *         by other task.
   */
  static bool TransportTracksTask(GeantTaskData *td, EventSet *workload);
  
  static
  FeederResult PreloadTracksForStep(GeantTaskData *td);
  
  static
  int SteppingLoop(GeantTaskData *td, bool flush);

  static
  int FlushOneLane(GeantTaskData *td);

  /**
   * @brief Function that provides transport tracks in coprocessor
   *
   * @param arg Arguments to be passed in the function
   */
  static void *TransportTracksCoprocessor(GeantPropagator *prop,TaskBroker *broker);

  /** @brief Function that provides waiting of workers */
  void WaitWorkers();
  
  int ShareBaskets(WorkloadManager *other);

private:
  /**
   * @brief Copy constructor for WorkloadManager
   * @details Still not implemented
   */
  WorkloadManager(const WorkloadManager &);

  /**
   * @brief Operator =
   * @todo Still not implemented
   */
  WorkloadManager &operator=(const WorkloadManager &);
};

} // GEANT_IMPL_NAMESPACE
} // Geant
#endif
