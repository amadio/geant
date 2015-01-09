//===--- GeantRqueue.h - Geant-V --------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantRqueue.h
 * @brief Definition of general queue for runnables.
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_RQUEUE
#define GEANT_RQUEUE

#include <deque>
#include "TCondition.h"
#include "TMutex.h"

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef GEANT_VRUNNABLE
#include "GeantVRunnable.h"
#endif

/** @brief GeantRqueue - A general queue for runnables. */
class GeantRqueue : public TObject {
private:
  std::deque<GeantVRunnable *> fQueue; /** Queue of runnables */
  mutable TMutex fMutex;               /** Mutex for the queue */
  TCondition fCondition;               /** Condition to take next runnable */
  Int_t fMax;                          /** Maximum allowed size */
  Int_t fN;                            /** Numer of objects in the queue */
  Int_t fNp;                           /** Number of prioritized objects */

public:

  /**
   * @brief Push function
   * @details Blocking push. The caller is quarantined if the queue is full
   * 
   * @param data GeantVRunnable data object
   * @param priority Priority parameter (by default false)
   */
  void Push(GeantVRunnable *data, Bool_t priority = kFALSE);
  
  /**
   * @brief Function that check possibility of push
   * @details Try push. If the queue is full return false and the data object is not injected
   * 
   * @param data GeantVRunnable data object
   * @param priority Priority parameter (by default false)
   * 
   */
  Bool_t TryPush(GeantVRunnable *data, Bool_t priority = kFALSE);
  
  /**
   * @brief Function that run next runnable
   * @details Run next runnable in the queue, possibly. Blocking if the queue is empty.
   */
  GeantVRunnable *RunNext();
  
  /**
   * @brief Function that try to run next runnable
   * @details Try to run next runnable in the queue. Non-blocking even if the queue is empty
   */
  GeantVRunnable *TryRunNext();

// Extract objects from the queue
