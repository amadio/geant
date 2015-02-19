//===--- GeantTaskManager.h - Geant-V ---------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTaskManager.h
 * @brief Implementation of a fixed-size block of user hits stored contiguously in memory
 * The class is used to efficiently manage hit user data in a multithreaded environment.
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TASKMANAGER
#define GEANT_TASKMANAGER

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

#include <vector>
#include "sync_objects.h"

/**
 * @brief Function that execute task function
 * 
 * @param data Data to be executed
 */
typedef void (*ExecuteTaskFunc_t)(void *data);

class TThread;

/**
 * @brief Class GeantTaskManager
 * @details A general concurrent task manager steered by data availability.
 * Use: RegisterFunction() method to register an arbitrary number of functions of type:
 * @code void Function(void *data)
 * @endcode that can process different data types. One has to use the
 * value returned by the RegisterFunction to identify the
 * work type when providing new workload to the task manager.
 */
class GeantTaskManager : public TNamed {
private: 
  Int_t fNactive;                      /** Number of active workers */
  Bool_t fToClean;                     /** Flag for requested cleanup */
  vector<TThread *> fThreads;          /** vector of threads */
  vector<ExecuteTaskFunc_t> exec_task; /** Vector of tasks of type void Function(void *data) */
  tcqueue fQueue;                      /** The work queue */
  Bool_t HasTask(ExecuteTaskFunc_t task) const;

public:

  /**
   * @brief GeantTaskManager constructor 
   * 
   * @param name Name for task manager
   */
  GeantTaskManager(const char *name);

  /** @brief Destructor GeantTaskManager */
  virtual ~GeantTaskManager() {}

  /**
   * @brief Function that add workers
   * 
   * @param nworkers Quantity of workers that should be added
   */
  Int_t AddWorkers(Int_t nworkers);

  /**
   * @brief Function that delete workers
   * 
   * @param nworkers Quantity of workers that should be deleted
   * @param sync Parameter for next synchronization (by default true)
   */
  void RemoveWorkers(Int_t nworkers, Bool_t sync = true);

  /**
   * @brief Function that stop workers
   * 
   * @param sync Parameter for next synchronization (by default true)
   */
  void StopWorkers(Bool_t sync = true);

  /**
   * @brief Function that return number tasks
   * @return Number of executed tasks
   */
  Int_t GetNtasks() const { return exec_task.size(); }

  /**
   * @brief Function that return size of vector of threads
   * @return Size of vector of threads
   */
  Int_t GetNthreads() const { return fThreads.size() }

  /**
   * @brief Function that return number of active workers
   * @return Number of active workers
   */
  Int_t GetNactive() const { return fNworkers; }

  /**
   * @brief Process loop function
   * 
   * @param arg Argument for processing
   */
  static void ProcessLoop(void *arg);

  /**
   * @brief Mark function for removal
   * 
   * @param arg Argument for removal
   */
  static void MarkForRemoval(void *arg);

  /**
   * @brief Function that register tasks
   * @details Main user methods. Register new task. Return task index.
   * 
   * @param task Task to register
   */
  Int_t RegisterTask(ExecuteTaskFunc_t task);
  
  /**
   * @brief Function that provides addition data to processing queue
   * @details Add new data chunk to the processing queue. The task index to process the
   * data must match the return value by RegisterTask
   * 
   * @param data Data that should be processed
   * @param task_id Task ID
   */
  void ToProcess(void *data, Int_t task_id);

  ClassDef(GeantTaskManager, 0) // A concurrent task manager
};
#endif
