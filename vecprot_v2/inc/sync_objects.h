//===--- sync_objects.h - Geant-V -------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file sync_objects.h
 * @brief Definition of synchronisation of objects in Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifdef USE_ROOT
#ifndef GEANT_SYNCOBJECTS
#define GEANT_SYNCOBJECTS
#include <deque>
#include <cassert>
#if __cplusplus >= 201103L
#include <atomic>
#endif
#include "TCondition.h"
#include "TMutex.h"

using namespace std;
class TStopwatch;
class TObject;

/** @struct TimeCounter */
struct TimeCounter {
  int nthreads; /** Number of threads */
  TStopwatch *timer;
  double stamp;
  double realtime[100];
  
  /**
   * @brief TimeCounter parameterized constructor
   * 
   * @param measure_time Measurement of time
   */
  TimeCounter(bool measure_time);

  /**
   * @brief TimeCounter copy constructor
   * @todo Not implemented
   */
  TimeCounter(const TimeCounter &);
  
  /** @brief TimeCounter destructor */
  ~TimeCounter();

  /** @brief Operator = */
  TimeCounter &operator=(const TimeCounter &); 

  /** @brief Operator ++ */
  TimeCounter &operator++();

  /** @brief Operator -- */
  TimeCounter &operator--();

  /** @brief Print function */
  void Print();
};

/** @brief Class for concurrent queue */
class concurrent_queue {
private:
  deque<TObject *> the_queue;
  mutable TMutex the_mutex;
  TCondition the_condition_variable;
  TimeCounter *the_counter;
  int nobjects;
  int npriority;
  
  /**
   * @brief Copy constructor for concurrent queue
   * @todo  Still not implemented
   */
  concurrent_queue(const concurrent_queue &); 
  /**
   * @brief Operator =
   * @todo  Still not implemented
   */
  concurrent_queue &operator=(const concurrent_queue &); // not implemented
public:
  
  /**
   * @brief Concurrent queue constructor
   * 
   * @param counter Counter for queue (by default false)
   */
  concurrent_queue(bool counter = false);
  
  /** @brief Concurrent queue destructor */
  ~concurrent_queue();
  
  /** @brief Function that assigned workers */
  int assigned_workers() const { return the_counter->nthreads; }

  /**
   * @brief Push function
   * 
   * @param data Data to be pushed
   * @param priority Priority of events (by default false)
   */
  void push(TObject *data, bool priority = false);

  /** @brief Function for sizing */
  int size() const;

  /** @brief Asynchronous function for sizing */
  int size_async() const { return nobjects; }

  /** @brief Function for empty */
  bool empty() const;

  /** @brief Asynchronized function for empty */
  bool empty_async() const { return nobjects == 0; }

  /** @brief Function for wait and pop */
  TObject *wait_and_pop();

  /**
   * @brief Function for maximization of wait and pop
   * 
   * @param nmax Maximum number
   * @param n Number of elements
   * @param array Array function
   */
  TObject *wait_and_pop_max(unsigned int nmax, unsigned int &n, TObject **array);

  /**
   * @brief Function for pop of many objects
   * 
   * @param n Number of objects
   * @param array Array of objects
   */
  void pop_many(unsigned int n, TObject **array);

  /** @brief Size of priority objects in queue */
  int size_priority() const { return npriority; }

  /** @brief Size of objects in queue */
  int size_objects() const { return nobjects; }

  /** @brief Print function */
  void Print();
};
#endif
/** @brief Reference counted atomic pointer */
template <class T> class ref_ptr {
public:

  T *fObjPtr; /** Object pointer */
#if __cplusplus >= 201103L
  std::atomic_flag fAcqLock; /** Memory barrier to lock acquiring */
  std::atomic_flag fXcgLock; /** Memory barrier to lock exchanging */
  std::atomic_int fRef;      /** Reference counter */
#endif
public:

  /** @brief Reference counted atomic pointer constructor */
  ref_ptr() : fObjPtr(0), fAcqLock(false), fXcgLock(false), fRef(0) {}

  /**
   * @brief Reference counted atomic pointer parametrized constructor
   * 
   * @param objptr Object pointer
   */
  ref_ptr(T *objptr) : fObjPtr(objptr), fAcqLock(false), fXcgLock(false), fRef(0) {}
  
  /**
   * @brief Set reference counted atomic pointer function
   * 
   * @param obj Object to which reference counted atomic pointer should be set
   */
  void set(T *obj) { fObjPtr = obj; }

  /** @brief Acquire function */
  T *acquire();

  /** @brief Clear function for memory barrier to lock acquiring */
  void clear_acq() { fAcqLock.clear(); }

  /** @brief Clear function for memory barrier to lock exchanging */
  void clear_xcg() { fXcgLock.clear(); }

  /** @brief Release function */
  void release();

  /** @brief Release and wait function */
  void release_and_wait();

  /**
   * @brief Replace and wait function
   * 
   * @param expected Expected pointer
   * @param newptr New pointer
   */
  bool replace_and_wait(T *expected, T *newptr);
};

/**
 * @brief Acquire function
 * @details Acquire the pointer while issuing a memory barrier for other possible clients
 */
template <class T> inline T *ref_ptr<T>::acquire() {
  while (fAcqLock.test_and_set()) {
  };                // barrier here for other threads
                    //__mutexed code start
  T *cpy = fObjPtr; // private snapshot
  fRef++;
  //__mutexed code end
  fAcqLock.clear();
  return cpy;
}

/**
 * @brief Release function
 * @details Release the copy of the pointer
 */
template <class T> inline void ref_ptr<T>::release() {
  fRef--;
  assert(fRef.load() >= 0);
}

/**
 * @brief Release and wait function
 * @details Release the copy of the pointer, but wait for everybody else to do it.
 * This acts as a spinlock.
 */
template <class T> inline void ref_ptr<T>::release_and_wait() {
  fRef--;
  while (fRef.load() > 0) {
  };
  assert(fRef.load() >= 0);
}

/**
 * @brief Replace and wait function
 * @details Atomically replace the pointer only if the old equals expected 
 * and only when it is not used anymore. Returns true if the replacement 
 * succeeded, in which case acquiring is blocked. The lock must be released
 * by the calling thread using clear() afterwards
 * 
 * @param expected Expected pointer
 * @param newptr New pointer
 */
template <class T> inline bool ref_ptr<T>::replace_and_wait(T *expected, T *newptr) {
  if (fXcgLock.test_and_set()) {
    // someone already passed for the same object
    release();
    while (fXcgLock.test_and_set()) {
    };
    fXcgLock.clear();
    assert(fRef.load() >= 0);
    return false;
  }
  // If someone has stolen the pointer, exit
  if (fObjPtr != expected) {
    // you thief!
    // Won't rely on expected pointer again
    release();
    assert(fRef.load() >= 0);
    fXcgLock.clear();
    return false;
  }
  // Thou shalt not pass
  fObjPtr = newptr;
  // go on, there's nothing to steal anymore
  // Lock acquiring new references
  while (fAcqLock.test_and_set()) {
  };
  release_and_wait();
  // now you're all mine...
  fAcqLock.clear();
  fXcgLock.clear();
  return true;
}

/** 
 * @brief Class basepipe
 * @details It has a dequeue<T> plus a function Process(T*) provided by the user
 */
class basepipe {
protected:
  int nobjects;  /** Number of objects in the queue */
  int npriority; /** Number of prioritized objects */
  int priority;  /** Assigned pipe priority */

public:

  /** @brief Basepipe constructor */
  basepipe() : nobjects(0), npriority(0), priority(0) {}

  /**
   * @brief basepipe destructor
   */
  virtual ~basepipe() {}
  
  /**
   * @brief Function that empty basepipe
   */
  bool empty() const { return (nobjects == 0); }

  /**
   * @brief Function that provides number of priority objects in basepipe
   */
  int size_priority() const { return npriority; }

  /**
   * @brief Function that provides number of objects in basepipe
   */
  int size() const { return nobjects; }

  /**
   * @brief Virtual function that provides process
   */
  virtual void process() = 0;
};

/**
 * @brief Templated class that provides workpipes public basepipe
 */
template <class T> class workpipe : public basepipe {
  typedef void *(*ProcessFunc_t)(T *data);
  ProcessFunc_t the_function;        /** Fuction to process this data */
  deque<T *> the_queue;              /** Double-ended data queue */
  mutable TMutex the_mutex;          /** General mutex for the queue */
  TCondition the_condition_variable; /** Condition */
  
  /**
   * @brief Function for processing user function in case of any object
   *  available in queue for workpipe
   */
  bool process_if_any();

public:

  /**
   * @brief Workpipe constructor
   * 
   * @param func Process function
   */
  workpipe(ProcessFunc_t func)
      : basepipe(), the_function(func), the_queue(), the_mutex(),
        the_condition_variable(&the_mutex) {}

  /** @brief Workpipe destructor */
  ~workpipe() {}

  /**
   * @brief Push function
   * 
   * @param data Data to be pushed
   * @param priority Priority state
   */
  void push(T *data, bool priority = false);

  /**
   * @brief Function for processing user function for workpipe
   */
  void process();
};

/**
 * @details Push an pointer of type T* in the queue. If pushed with priority, the pointer
 * is put at the bask of the queue, otherwise to the front.
 */
template <class T> void workpipe<T>::push(T *data, bool pr) {
  the_mutex.Lock();
  nobjects++;
  if (pr) {
    the_queue.push_back(data);
    npriority++;
  } else
    the_queue.push_front(data);
  the_condition_variable.Signal();
  the_mutex.UnLock();
}

/**
 * @details Gets the back object from the queue. Wait if none is available. Call the
 * user process function for the first available object.
 */
template <class T> void workpipe<T>::process() {
  the_mutex.Lock();
  while (the_queue.empty())
    the_condition_variable.Wait();
  T *popped_value = the_queue.back();
  the_queue.pop_back();
  nobjects--;
  if (npriority > 0)
    npriority--;
  the_mutex.UnLock();
  // Execute the user function
  the_function(popped_value);
}

/**
 * @details If any object is available in the queue, pop it and call the
   * user process function for it. Returns true if the processing was called
 */
template <class T> bool workpipe<T>::process_if_any() {
  the_mutex.Lock();
  if (!nobjects) {
    the_mutex.UnLock();
    return false;
  }
  T *popped_value = the_queue.back();
  the_queue.pop_back();
  nobjects--;
  if (npriority > 0)
    npriority--;
  the_mutex.UnLock();
  // Execute the user function
  the_function(popped_value);
  return true;
}

// class workflow

#endif
