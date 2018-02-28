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
 * @param m  Mutex
 * @param cv  Condition variable
 * @param start  Starting atomic variable
 */
struct condition_locker {
  std::mutex m;
  std::condition_variable cv;
  std::atomic<bool> start;

  condition_locker() : m(), cv(), start(false) {}

  /**
   * @brief  Wait function for controlled threads
   * @details Simple function based on idea that controlled threads have
   * to call the Wait() function bofore starting waiting thread(s)
   */
  void Wait()
  {
    std::unique_lock<std::mutex> lk(m);
    while (!start.load())
      cv.wait(lk);
    start.store(false);
    lk.unlock();
    cv.notify_one();
  };

  /**
   * @brief Function of starting one waiting thread
   * @details Controller starts one thread after calling by controlled thread function Wait()
   */
  void StartOne()
  {
    std::unique_lock<std::mutex> lk(m);
    start.store(true);
    cv.notify_one();
  };

  /**
   * @brief Function of starting of all waiting threads
   * @details Controller starts all threads after calling by controlled thread function Wait()
   */
  void StartAll()
  {
    std::unique_lock<std::mutex> lk(m);
    start.store(true);
    cv.notify_all();
  };
};
#endif
