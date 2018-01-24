//===--- mpmc_bounded_queue.h - Geant-V -------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file mpmc_bounded_queue.h
 * @brief Implementation of MPMC bounded queue for Geant-V prototype 
 * @author Andrei Gheata
 * based on http://www.1024cores.net/home/lock-free-algorithms/queues/bounded-mpmc-queue
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_MPMC_BOUNDED_QUEUE
#define GEANT_MPMC_BOUNDED_QUEUE
#include <atomic>
#include <cassert>

typedef std::size_t size_t;

/**
 * @brief Class MPMC bounded queue
 * 
 * @param buffer_size Buffer size for queue
 * @tparam T Type of objects
 */
template <typename T> class mpmc_bounded_queue {
public:

  /**
   * @brief MPMC bounded queue constructor
   * 
   * @param buffer_size Buffer size for queue
   */
  mpmc_bounded_queue(size_t buffer_size)
      : buffer_(new cell_t[buffer_size]), buffer_mask_(buffer_size - 1), enqueue_pos_(0),
        dequeue_pos_(0), nstored_(0) {
    assert((buffer_size >= 2) && ((buffer_size & (buffer_size - 1)) == 0));
    for (size_t i = 0; i != buffer_size; i += 1)
      buffer_[i].sequence_.store(i, std::memory_order_relaxed);
    enqueue_pos_.store(0, std::memory_order_relaxed);
    dequeue_pos_.store(0, std::memory_order_relaxed);
    nstored_.store(0, std::memory_order_relaxed);
  }

  /** @brief MPMC bounded queue destructor */
  ~mpmc_bounded_queue() { delete[] buffer_; }
  
  /** @brief Size function */
  inline
  size_t size() const { return nstored_.load(std::memory_order_relaxed); }
  
  /**
   * @brief MPMC enqueue function
   * 
   * @param data Data to be enqueued
   */
  inline
  bool enqueue(T const &data) {
    cell_t *cell;
    size_t pos = enqueue_pos_.load(std::memory_order_relaxed);
    for (;;) {
      cell = &buffer_[pos & buffer_mask_];
      size_t seq = cell->sequence_.load(std::memory_order_acquire);
      intptr_t dif = (intptr_t)seq - (intptr_t)pos;
      if (dif == 0) {
        if (enqueue_pos_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed))
          break;
      } else if (dif < 0)
        return false;
      else
        pos = enqueue_pos_.load(std::memory_order_relaxed);
    }
    cell->data_ = data;
    nstored_++;
    cell->sequence_.store(pos + 1, std::memory_order_release);
    return true;
  }
  
  /**
   * @brief MPMC dequeue function
   * 
   * @param data Data to be dequeued
   */
  inline
  bool dequeue(T &data) {
    cell_t *cell;
    size_t pos = dequeue_pos_.load(std::memory_order_relaxed);
    for (;;) {
      cell = &buffer_[pos & buffer_mask_];
      size_t seq = cell->sequence_.load(std::memory_order_acquire);
      intptr_t dif = (intptr_t)seq - (intptr_t)(pos + 1);
      if (dif == 0) {
        if (dequeue_pos_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed))
          break;
      } else if (dif < 0)
        return false;
      else
        pos = dequeue_pos_.load(std::memory_order_relaxed);
    }
    data = cell->data_;
    nstored_--;
    cell->sequence_.store(pos + buffer_mask_ + 1, std::memory_order_release);
    return true;
  }

  /**
   * @brief MPMC dequeue function trying to fetch into atomic variable, if this contains the expected value
   * 
   * @param data Data to be dequeued
   */
  inline
  bool dequeue(std::atomic<T> &data, T &expected) {
    cell_t *cell;
    size_t pos = dequeue_pos_.load(std::memory_order_relaxed);
    for (;;) {
      cell = &buffer_[pos & buffer_mask_];
      size_t seq = cell->sequence_.load(std::memory_order_acquire);
      intptr_t dif = (intptr_t)seq - (intptr_t)(pos + 1);
      if (dif == 0) {
        if (!data.compare_exchange_strong(expected, cell->data_))
          return false;
        if (dequeue_pos_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed))
          break;
      } else if (dif < 0)
        return false;
      else
        pos = dequeue_pos_.load(std::memory_order_relaxed);
    }
    expected = cell->data_;
    nstored_--;
    cell->sequence_.store(pos + buffer_mask_ + 1, std::memory_order_release);
    return true;
  }


private:

  /** @struct cell_t */
  struct cell_t {
    std::atomic<size_t> sequence_;
    T data_;

    /**
     * @brief Cell constructor
     * 
     * @param sequence_(0) Sequence 
     * @param data_() Data
     */
    cell_t() : sequence_(0), data_() {}
  };

  static size_t const cacheline_size = 64;
  typedef char cacheline_pad_t[cacheline_size];
  cacheline_pad_t pad0_;
  cell_t *const buffer_;
  size_t const buffer_mask_;
  cacheline_pad_t pad1_;
  std::atomic<size_t> enqueue_pos_;
  cacheline_pad_t pad2_;
  std::atomic<size_t> dequeue_pos_;
  cacheline_pad_t pad3_;
  std::atomic<size_t> nstored_;
  cacheline_pad_t pad4_;
  
  /** @brief MPMC bounded queue copy constructor */
  mpmc_bounded_queue(mpmc_bounded_queue const &);

  /** @brief Operator = */
  void operator=(mpmc_bounded_queue const &);
};
#endif
