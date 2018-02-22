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

#include "Propagator.h"

#include "TaskData.h"

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
  Propagator* fPropagator = nullptr; /** propagator,... */
  int fNthreads = 0;                      /** Number of managed threads */
  bool fStarted = false;                  /** Start flag */
  bool fStopped = false;                  /** Stop flag */
  priority_queue<Basket *> *fDoneQ = nullptr; /** Thread "all work done" queue */
  std::vector<std::thread> fListThreads;  /** Vector of threads */

  TaskBroker *fBroker = nullptr;          /** Pointer to the coprocessor broker, this could be made a collection. */

  /**
   * @brief WorkloadManager parameterized constructor
   *
   * @param  nthreads Number of threads foe workload manager
   */
  WorkloadManager(int nthreads, Propagator* prop) : fPropagator(prop), fNthreads(nthreads)
  {
    fDoneQ = new priority_queue<Basket *>(1 << 10, nthreads);
  }

public:
  /** @brief WorkloadManager destructor */
  virtual ~WorkloadManager() 
  {
    delete fDoneQ;
  }

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
  static WorkloadManager *NewInstance(Propagator *prop= nullptr, int nthreads = 0);

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
  TaskBroker *GetTaskBroker() { return fBroker; }

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
  
  static
  FeederResult PreloadTracksForStep(TaskData *td);
  
  static
  int SteppingLoop(TaskData *td, bool flush);

  static
  int FlushOneLane(TaskData *td);

  /** @brief Function that provides waiting of workers */
  void WaitWorkers();

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
