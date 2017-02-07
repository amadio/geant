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

#include <mutex>
#include <condition_variable>

namespace Geant {

/** @brief Class of priority queue */
template <typename T> class priority_queue {
public:

  /**
   * @brief Priority_queue constructor
   * 
   * @param buffer_size Buffer size for queue
   */
  priority_queue(size_t buffer_size, size_t threshold = 0)
      : n_waiting_(0), mutex_(), cv_(), main_q_(buffer_size),
        priority_q_(buffer_size), threshold_(threshold), counter_(0) {}
  
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
  inline bool push(T const &data, bool priority = false);
  
  /**
   * @brief Forced push function
   * 
   * @param data Data to be pushed
   * @param priority Priority state (by default false)
   */
  inline void push_force(T const &data, bool priority = false);
 
   /** @brief Size function of number of enqueued objects */
  inline size_t size_async() const;
  
  /** @brief Size function of priority objects in queue */
  inline size_t size_priority() const;
  
 /** @brief Function that returns numbere of enqueued objects in the priority queue */
  inline size_t size_objects() const;
  
  /** @brief Function that empty number of enqueued objects */
  inline bool empty_async() const;
  
  /** @brief Size function for number of enqueued objects*/
  inline int size() const { return size_async(); }

   /** @brief Size function for number of enqueued objects*/
  inline unsigned int status() const;
 
  /** @brief Empty function for number of enqueued objects */
  inline bool empty() const { return empty_async(); }
  
  /**
   * @brief Pop function
   * 
   * @param data Data that should be popped
   */
  inline bool try_pop(T &data);
  
  /**
   * @brief Wait and pop function
   * 
   * @param data Data that should be popped
   */
  inline void wait_and_pop(T &data);
      
private:
  std::atomic_int n_waiting_;  /** number of threads waiting */
  std::mutex mutex_;           /** mutex for starvation mode */
  std::condition_variable cv_; /** condition variable */

  mpmc_bounded_queue<T> main_q_;
  mpmc_bounded_queue<T> priority_q_;

  size_t threshold_;           /** Threshold for the normal filling */
  std::atomic<size_t> counter_; /** Queue content counter */
};

/**
 * @details Push one element in the queue. Non-blocking, unless there is not enough
 * allocated space.
 * @return List of pushed objects
 */
template <typename T> 
bool priority_queue<T>::push(T const &data, bool priority) {
  if ((priority && priority_q_.enqueue(data)) || main_q_.enqueue(data)) {
    counter_.fetch_add(1, std::memory_order_relaxed);
    if (n_waiting_.load(std::memory_order_relaxed))
      cv_.notify_one();
    return true;
  }
  return false;
}

template <typename T> 
void priority_queue<T>::push_force(T const &data, bool priority) {
  while ( !push(data, priority) ) {}
}

/** @todo Add details */
template <typename T>
size_t priority_queue<T>::size_async() const {
  return ( counter_.load(std::memory_order_relaxed) );
}

/** @todo  Add details */
template <typename T>
size_t priority_queue<T>::size_priority() const {
  return priority_q_.size();
}

/**
 * @brief Function that returns number of enqueued objects in the priority queue
 */
template <typename T>
size_t priority_queue<T>::size_objects() const {
  return main_q_.size();
}

/** @brief  Asynchronous check if queue is empty */
template <typename T> 
bool priority_queue<T>::empty_async() const {
  return ( counter_.load(std::memory_order_relaxed) == 0 );
}

/** @brief Queue status getter.
  *  @return Queue status: 0=empty, 1=depleted, 2=normal
  */
template <typename T> 
unsigned int priority_queue<T>::status() const { 
  size_t counter = counter_.load(std::memory_order_relaxed);
  if (counter == 0) return 0;
  if (counter < threshold_) return 1;
  return 2;
}

/** 
 * @brief Unblocking retrieval of one element
 * @return  Success of operation
 */
template <typename T>
bool priority_queue<T>::try_pop(T &data) {
  if ( priority_q_.dequeue(data) || main_q_.dequeue(data) ) {
    counter_.fetch_sub(1, std::memory_order_relaxed);
    return true;
  }
  return false;
}

/** @brief Blocking pop operation retrieving one element from the queue */
template <typename T>
void priority_queue<T>::wait_and_pop(T &data) {
  if (!n_waiting_.load() && try_pop(data))
    return;
  // We have to put the thread to sleep
  n_waiting_++;
  std::unique_lock<std::mutex> lk(mutex_);
  while (size_async() == 0)
    cv_.wait(lk);
  while (!try_pop(data)) {}
  n_waiting_--;
  lk.unlock();
}
} // namespace Geant

#endif
