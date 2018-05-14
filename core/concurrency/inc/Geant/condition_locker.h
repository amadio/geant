//===--- condition_locker.h - Geant-V ---------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file condition_locker.h
 * @brief Implementation of condition locker
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_CONDITIONLOCKER
#define GEANT_CONDITIONLOCKER

#include <mutex>
#include <condition_variable>

/**
 * @struct condition_locker
 * @brief A simple semaphore based on condition.
 * @details The controlled threads have to call the Wait() function and the controller
 *  can wake up one waiting thread at a time calling StartOne() method, or all by calling StartAll().
 *
 * @param fMutex Mutex
 * @param fCondition Condition variable
 * @param fCanStart Starting atomic variable
 */
struct condition_locker {
  std::atomic_int fNwait;             ///< number of threads waiting on the semaphore
  std::mutex fMutex;                  ///< semaphore mutex
  std::condition_variable fCondition; ///< semaphore condition
  std::atomic<bool> fCanStart;        ///< starter atomic condition

  condition_locker() : fNwait(0), fMutex(), fCondition(), fCanStart(false) {}

  /** @brief Getter for number of threads waiting */
  int GetNwait() const { return fNwait.load(); }

  /**
   * @brief  Wait function for controlled threads
   * @details Simple function based on idea that controlled threads have
   * to call the Wait() function bofore starting waiting thread(s)
   */
  void Wait()
  {
    std::unique_lock<std::mutex> lk(fMutex);
    fNwait++;
    while (!fCanStart.load())
      fCondition.wait(lk);
    fCanStart.store(false);
    fNwait--;
    lk.unlock();
    // fCondition.notify_one();
  };

  /**
   * @brief Function of starting one waiting thread
   * @details Controller starts one thread after calling by controlled thread function Wait()
   */
  void StartOne()
  {
    std::unique_lock<std::mutex> lk(fMutex);
    fCanStart.store(true);
    fCondition.notify_one();
  };

  /**
   * @brief Function of starting of all waiting threads
   * @details Controller starts all threads after calling by controlled thread function Wait()
   */
  void StartAll()
  {
    std::unique_lock<std::mutex> lk(fMutex);
    fCanStart.store(true);
    fCondition.notify_all();
  };
};
#endif
