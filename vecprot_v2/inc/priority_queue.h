//===--- priority_queue.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file priority_queue.h
 * @brief Implementation of priority queue used in Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PRIORITY_QUEUE
#define GEANT_PRIORITY_QUEUE

#ifndef GEANT_MPMC_BOUNDED_QUEUE
#include "mpmc_bounded_queue.h"
#endif

#if __cplusplus >= 201103L
#include <mutex>
#include <condition_variable>
#endif

namespace Geant {

/** @brief Class of priority queue */
template <typename T> class priority_queue {
public:

  /**
   * @brief Priority_queue constructor
   * 
   * @param buffer_size Buffer size for queue
   */
  priority_queue(size_t buffer_size)
      : n_waiting_(0), countdown_(0), mutex_(), cv_(), main_q_(buffer_size),
        priority_q_(buffer_size) {}
  
  /**
   * @brief  Priority_queue destructor
   */
  ~priority_queue() {}
  
  /**
   * @brief Push function
   * 
   * @param data Data to be pushed
   * @param priority Priority state (by default false)
   */
  bool push(T const &data, bool priority = false);
  
  /** @brief Size function of number of enqueued objects */
  size_t size_async() const;
  
  /** @brief Size function of priority objects in queue */
  size_t size_priority() const;
  
 /** @brief Function that returns numbere of enqueued objects in the priority queue */
  size_t size_objects() const;
  
  /** @brief Function that empty number of enqueued objects */
  bool empty_async() const;
  
  /** @brief Size function for number of enqueued objects*/
  int size() const { return size_async(); }
  
  /** @brief Empty function for number of enqueued objects */
  bool empty() const { return empty_async(); }
  
  /**
   * @brief Pop function
   * 
   * @param data Data that should be popped
   */
  bool try_pop(T &data);
  
  /**
   * @brief Wait and pop function
   * 
   * @param data Data that should be popped
   */
  void wait_and_pop(T &data);
  
  /**
   * @brief Funtion that get countdown value
   */
  size_t get_countdown() const { return countdown_.load(std::memory_order_relaxed); }
  
  /**
   * @brief Funtion that reset countdown value
   */
  void reset_countdown() { countdown_.store(-1, std::memory_order_relaxed); }
  
  /**
   * @brief Funtion that set countdown value
   * 
   * @param n Countdown counter for extracted objects
   */
  void set_countdown(size_t n) { countdown_.store(n, std::memory_order_relaxed); }
  
  /**
   * @brief ??????
   */
  size_t n_ops() const { return 0; }

private:
#if __cplusplus >= 201103L
  std::atomic_int n_waiting_;  /** number of threads waiting */
  std::atomic_int countdown_;  /** Countdown counter for extracted objects */
  std::mutex mutex_;           /** mutex for starvation mode */
  std::condition_variable cv_; /** condition variable */
#endif
  mpmc_bounded_queue<T> main_q_;
  mpmc_bounded_queue<T> priority_q_;
};

/**
 * @details Push one element in the queue. Non-blocking, unless there is not enough
 * allocated space.
 * @return List of pushed objects
 */
template <typename T> bool priority_queue<T>::push(T const &data, bool priority) {
  bool pushed = false;
  if (priority)
    pushed = priority_q_.enqueue(data);
  if (!pushed)
    pushed = main_q_.enqueue(data);
  if (n_waiting_.load(std::memory_order_relaxed))
    cv_.notify_one();
  return pushed;
}

/** @todo Add details */
template <typename T> size_t priority_queue<T>::size_async() const {
  return (main_q_.size() + priority_q_.size());
}

/** @todo  Add details */
template <typename T> size_t priority_queue<T>::size_priority() const {
  return priority_q_.size();
}

/**
 * @brief Function that returns number of enqueued objects in the priority queue
 */
template <typename T> size_t priority_queue<T>::size_objects() const {
  return main_q_.size();
}

/** @todo  Add details */
template <typename T> bool priority_queue<T>::empty_async() const {
  return (size_async() == 0);
}

/** 
 * @todo  Add details
 * @return  List of popped elements
 */
template <typename T> bool priority_queue<T>::try_pop(T &data) {
  bool popped = priority_q_.dequeue(data);
  if (!popped)
    popped = main_q_.dequeue(data);
  if (popped && countdown_.load() > 0)
    countdown_--;
  return popped;
}

/** @todo  Add details */
template <typename T> void priority_queue<T>::wait_and_pop(T &data) {
  if (!n_waiting_.load() && try_pop(data))
    return;
  // We have to put the thread to sleep
  n_waiting_++;
  std::unique_lock<std::mutex> lk(mutex_);
  while (size_async() == 0)
    cv_.wait(lk);
  while (!try_pop(data)) {
  }
  n_waiting_--;
  if (countdown_.load() > 0)
    countdown_--;
  lk.unlock();
}
} // namespace Geant

#endif
