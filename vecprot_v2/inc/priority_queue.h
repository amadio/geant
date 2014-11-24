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

template<typename T>
class priority_queue
{
public:
  priority_queue(size_t buffer_size) 
  : n_waiting_(0),
    countdown_(0),
    mutex_(),
    cv_(),
    main_q_(buffer_size),
    priority_q_(buffer_size) {}
  ~priority_queue()          {}

   bool               push(T const &data, bool priority=false);
   size_t             size_async() const;
   size_t             size_priority() const;
   size_t             size_objects() const;
   bool               empty_async() const;
   int                size() const {return size_async();}
   bool               empty() const {return empty_async();}
   bool               try_pop(T& data);
   void               wait_and_pop(T& data);
   size_t             get_countdown() const {return countdown_.load(std::memory_order_relaxed);}   
   void               reset_countdown() {countdown_.store(-1,std::memory_order_relaxed);}
   void               set_countdown(size_t n) {countdown_.store(n,std::memory_order_relaxed);}
   size_t             n_ops() const {return 0;}

private:
#if __cplusplus >= 201103L
  std::atomic_int         n_waiting_;   // number of threads waiting  
  std::atomic_int         countdown_;   // Countdown counter for extracted objects
  std::mutex              mutex_;       // mutex for starvation mode
  std::condition_variable cv_;          // condition variable
#endif  
  mpmc_bounded_queue<T>   main_q_;  
  mpmc_bounded_queue<T>   priority_q_;
};

template<typename T>
bool priority_queue<T>::push(T const &data, bool priority)
{
// Push one element in the queue. Non-blocking, unless there is not enough
// allocated space.
  bool pushed = false;
  if (priority) pushed = priority_q_.enqueue(data);
  if (!pushed)  pushed = main_q_.enqueue(data);
  if (n_waiting_.load(std::memory_order_relaxed)) cv_.notify_one();
  return pushed;
}

template<typename T>
size_t priority_queue<T>::size_async() const
{
// Returns numbere of enqueued objects
   return (main_q_.size()+priority_q_.size());
}

template<typename T>
size_t priority_queue<T>::size_priority() const
{
// Returns numbere of enqueued objects in the priority queue
   return priority_q_.size();
}

template<typename T>
size_t priority_queue<T>::size_objects() const
{
// Returns numbere of enqueued objects in the priority queue
   return main_q_.size();
}

template<typename T>
bool priority_queue<T>::empty_async() const
{
// Checks if the queue is empty
  return (size_async() == 0);
}

template<typename T>
bool priority_queue<T>::try_pop(T& data)
{
   bool popped = priority_q_.dequeue(data);
   if (!popped) popped = main_q_.dequeue(data);
   if (popped && countdown_.load()>0) countdown_--;
   return popped;
}

template<typename T>
void priority_queue<T>::wait_and_pop(T& data)
{
  if (!n_waiting_.load() && try_pop(data)) return;
  // We have to put the thread to sleep
  n_waiting_++;
  std::unique_lock<std::mutex> lk(mutex_);
  while (size_async()==0) cv_.wait(lk);
  while (!try_pop(data)) {}
  n_waiting_--;
  if (countdown_.load()>0) countdown_--;
  lk.unlock();
}

}

#endif
