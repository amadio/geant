#if __cplusplus >= 201103L
#include <atomic>
#endif
#include <cassert>

typedef std::size_t size_t;

template<typename T>
class mpmc_bounded_queue
{
public:
  mpmc_bounded_queue(size_t buffer_size)
    : buffer_(new cell_t [buffer_size])
    , buffer_mask_(buffer_size - 1)
    , enqueue_pos_(0)
    , dequeue_pos_(0)
    , nstored_(0)
  {
    assert((buffer_size >= 2) &&
      ((buffer_size & (buffer_size - 1)) == 0));
    for (size_t i = 0; i != buffer_size; i += 1)
      buffer_[i].sequence_.store(i, std::memory_order_relaxed);
    enqueue_pos_.store(0, std::memory_order_relaxed);
    dequeue_pos_.store(0, std::memory_order_relaxed);
    nstored_.store(0, std::memory_order_relaxed);
  }

  ~mpmc_bounded_queue()
  {
    delete [] buffer_;
  }

  size_t size() const
  {
    return nstored_.load(std::memory_order_relaxed);
  }
  
  bool enqueue(T const& data)
  {
    cell_t* cell;
    size_t pos = enqueue_pos_.load(std::memory_order_relaxed);
    for (;;)
    {
      cell = &buffer_[pos & buffer_mask_];
      size_t seq = 
        cell->sequence_.load(std::memory_order_acquire);
      intptr_t dif = (intptr_t)seq - (intptr_t)pos;
      if (dif == 0)
      {
        if (enqueue_pos_.compare_exchange_weak
            (pos, pos + 1, std::memory_order_relaxed))
          break;
      }
      else if (dif < 0)
        return false;
      else
        pos = enqueue_pos_.load(std::memory_order_relaxed);
    }
    cell->data_ = data;
    nstored_++;
    cell->sequence_.store(pos + 1, std::memory_order_release);
    return true;
  }

  bool dequeue(T& data)
  {
    cell_t* cell;
    size_t pos = dequeue_pos_.load(std::memory_order_relaxed);
    for (;;)
    {
      cell = &buffer_[pos & buffer_mask_];
      size_t seq = 
        cell->sequence_.load(std::memory_order_acquire);
      intptr_t dif = (intptr_t)seq - (intptr_t)(pos + 1);
      if (dif == 0)
      {
        if (dequeue_pos_.compare_exchange_weak
            (pos, pos + 1, std::memory_order_relaxed))
          break;
      }
      else if (dif < 0)
        return false;
      else
        pos = dequeue_pos_.load(std::memory_order_relaxed);
    }
    data = cell->data_;
    nstored_--;
    cell->sequence_.store
      (pos + buffer_mask_ + 1, std::memory_order_release);
    return true;
  }

private:
  struct cell_t
  {
#if __cplusplus >= 201103L
    std::atomic<size_t>   sequence_;
#endif
    T                     data_;
    cell_t(): sequence_(0), data_() {}
  };

  static size_t const     cacheline_size = 64;
#if __cplusplus >= 201103L
  typedef char            cacheline_pad_t [cacheline_size];
#else
  typedef char            cacheline_pad_t [64];
#endif
  cacheline_pad_t         pad0_;
  cell_t* const           buffer_;
  size_t const            buffer_mask_;
  cacheline_pad_t         pad1_;
#if __cplusplus >= 201103L
  std::atomic<size_t>     enqueue_pos_;
#endif
  cacheline_pad_t         pad2_;
#if __cplusplus >= 201103L
  std::atomic<size_t>     dequeue_pos_;
#endif
  cacheline_pad_t         pad3_;
#if __cplusplus >= 201103L
  std::atomic<size_t>     nstored_;
#endif
  cacheline_pad_t         pad4_;

  mpmc_bounded_queue(mpmc_bounded_queue const&);
  void operator = (mpmc_bounded_queue const&);
}; 
 
