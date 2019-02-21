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
#include "Geant/priority_queue.h"
#include "Geant/condition_locker.h"
#include "Geant/dcqueue.h"
#include "Geant/condition_locker.h"

#include "Geant/Propagator.h"

#include "Geant/TaskData.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantBasketMgr;
class Basket;
class GeantScheduler;
class TaskBroker;
class EventSet;
class RunManager;

/**
 * @brief WorkloadManager class
 * @details Main work manager class. This creates and manages all the worker threads,
 * has pointers to the synchronization objects, but also to the currently
 * transported baskets.
 */
class WorkloadManager {
protected:
  Propagator *fPropagator          = nullptr; /** propagator,... */
  int fNthreads                    = 0;       /** Number of managed threads */
  bool fStarted                    = false;   /** Start flag */
  bool fStopped                    = false;   /** Stop flag */
  priority_queue<Basket *> *fDoneQ = nullptr; /** Thread "all work done" queue */
  std::vector<std::thread> fListThreads;      /** Vector of threads */
  condition_locker fSemaphore;                /** Semaphore for starving threads */

  TaskBroker *fBroker = nullptr; /** Pointer to the coprocessor broker, this could be made a collection. */

  /**
   * @brief WorkloadManager parameterized constructor
   *
   * @param  nthreads Number of threads foe workload manager
   */
  WorkloadManager(int nthreads, Propagator *prop) : fPropagator(prop), fNthreads(nthreads)
  {
    fDoneQ = new priority_queue<Basket *>(1 << 10, nthreads);
  }

public:
  /** @brief WorkloadManager destructor */
  virtual ~WorkloadManager() { delete fDoneQ; }

  enum class FeederResult : char { kNone, kError, kWork, kStop };

  /** @brief Call Feeder (if needed) and check exit condition. */
  FeederResult CheckFeederAndExit();

  /** @brief Function that returns number of managed threads */
  GEANT_FORCE_INLINE
  int GetNthreads() const { return fNthreads; }

  /**
   * @brief Function that provides thread's "all work done" queue
   */
  priority_queue<Basket *> *DoneQueue() const { return fDoneQ; }

  /**
   * @brief Function that create workload manager instance
   *
   * @param nthreads Number of threads (by default 0)
   */
  static WorkloadManager *NewInstance(Propagator *prop = nullptr, int nthreads = 0);

  /** @brief Function that check stop flag */
  GEANT_FORCE_INLINE
  bool IsStopped() const { return fStopped; }

  /** @brief Function that provide stop process by setting Stop flag = True */
  GEANT_FORCE_INLINE
  void Stop() { fStopped = true; }

  /** @brief Stop all transport threads */
  GEANT_FORCE_INLINE
  void StopTransportThreads() { fStopped = true; }

  /** @brief  Setter for task broker */
  void SetTaskBroker(TaskBroker *broker);

  /** @brief  Setter for task broker */
  GEANT_FORCE_INLINE
  TaskBroker *GetTaskBroker() const { return fBroker; }

  /** @brief Function that initializes the threads/tasks used by the system */
  bool StartTasks();

  /** @brief Joins all threads at the end of processing */
  void JoinThreads();

#ifdef USE_ROOT
  /** @brief Function that starts ROOT application */
  static void *StartROOTApplication();
#endif

  /**
   * @brief Function that provides transporting tracks
   *
   * @param arg Arguments to be passed in the function
   */
  static void TransportTracksV3(Propagator *prop);

  /**
   * @brief Function that provides transporting tracks
   *
   * @param td Task data
   * @param workload Event set to be transported
   * @return True if workload completed. If false, the work will be completed
   *         by other task.
   */
  static bool TransportTracksTask(EventSet *workload, TaskData *td);

  /** This transport entry point emulates single track transport throughout stages */
  static void TransportTracksSingle(Propagator *prop);

  static FeederResult PreloadTracksForStep(TaskData *td);

  /** This is the standard V3 stepping loop. */
  static int SteppingLoop(TaskData *td, bool flush);

  /** This is the single track stepping loop */
  static int SteppingLoopSingle(TaskData *td);

  static int FlushOneLane(TaskData *td, bool share);

  /** @brief Function that provides waiting of workers */
  void WaitWorkers();

  /** @brief Get number of suspended workers */
  GEANT_FORCE_INLINE
  int GetNsuspended() const { return fSemaphore.GetNwait(); }

  /** @brief Suspend one worker */
  GEANT_FORCE_INLINE
  void Wait() { fSemaphore.Wait(); }

  /** @brief Start one suspended worker */
  GEANT_FORCE_INLINE
  void StartOne() { fSemaphore.StartOne(); }

  /** @brief Start all suspended workers */
  GEANT_FORCE_INLINE
  void StartAll() { fSemaphore.StartAll(); }

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

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
#endif
