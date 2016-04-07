//===--- Basketizer.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Basketizer.h
 * @brief Implementation of basketizer of T* pointers for Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_BASKETIZER
#define GEANT_BASKETIZER

#include <atomic>
#include <vector>
#include <cassert>
#include <iostream>
#include "Geant/Config.h"

/**
 * @brief Class Basketizer
 *
 * @param buffer_size Buffer size for queue
 * @tparam T Type of objects
 */
namespace Geant {

template <typename T> class BasketCounter {
  using size_t = std::size_t;

public:
  size_t fBsize;                // Basket size
  std::atomic<size_t> fNbooked; // Counter for booked slots
  std::atomic<size_t> fNfilled; // Counter for filled slots
  std::atomic_flag fLock;       // Lock for the basket counter
  BasketCounter() : fBsize(0), fNbooked(0), fNfilled(0), fLock() {
    // Constructor
    fLock.clear();
  }

  BasketCounter(short bsize) : fBsize(bsize), fNbooked(0), fNfilled(0), fLock() {
    // Constructor
    fLock.clear();
  }

  GEANT_INLINE
  void SetBsize(size_t bsize) { fBsize = bsize; }

  GEANT_INLINE
  size_t Bsize() { return (fBsize); }

  GEANT_INLINE
  size_t BookSlot() { return (fNbooked.fetch_add(1) + 1); }

  GEANT_INLINE
  size_t FillSlot(size_t &ibooked, T *address, T const data) {
    if (ibooked > fBsize) {
      // Block slot filling until basket is released
      while (Nbooked() > fBsize)
        assert(Nbooked() < 2*fBsize);
        ;
      ibooked -= fBsize;
    }
    *address = data;
    return (fNfilled.fetch_add(1) + 1);
  }

  GEANT_INLINE
  void ReleaseBasket() {
    fNfilled -= fBsize;
    fNbooked -= fBsize;
  }

  GEANT_INLINE
  void Clear() {
    fNfilled.store(0);
    fNbooked.store(0);
  }

  GEANT_INLINE
  size_t Nbooked() const { return (fNbooked.load()); }

  GEANT_INLINE
  size_t Nfilled() const { return (fNfilled.load()); }

  GEANT_INLINE
  void Lock() {
    while (fLock.test_and_set(std::memory_order_acquire))
      ;
  }

  GEANT_INLINE
  void Unlock() { fLock.clear(std::memory_order_release); }
};

template <typename T> class Basketizer {
public:
  using size_t = std::size_t;
  using BasketCounter_t = BasketCounter<T *>;

  /**
   * @brief Basketizer dummy constructor
   */
  Basketizer()
      : fBsize(0), fBmask(0), fLock(), fBuffer(0), fCounters(0), fBufferMask(0), fIbook(0), fNstored(0), fNbaskets(0) {
    fLock.clear();
  }
  /**
   * @brief Basketizer default constructor
   *
   * @param buffer_size Circular buffer size
   * @param basket_size Size of produced baskets
   */
  Basketizer(size_t buffer_size, unsigned int basket_size)
      : fBsize(0), fBmask(0), fLock(), fBuffer(0), fCounters(0), fBufferMask(0), fIbook(0), fNstored(0), fNbaskets(0) {
    fLock.clear();
    Init(buffer_size, basket_size);
  }

  /** @brief Basketizer destructor */
  ~Basketizer() {
    delete[] fBuffer;
    delete[] fCounters;
  }

  //____________________________________________________________________________
  /** @brief Initialize basketizer */
  void Init(size_t buffer_size, unsigned int basket_size) {
    // Make sure the requested size is a power of 2
    assert((buffer_size >= 2) && ((buffer_size & (buffer_size - 1)) == 0));
    assert((basket_size >= 2) && ((basket_size & (basket_size - 1)) == 0));
    fBsize = basket_size;
    fBuffer = new T *[buffer_size];
    fCounters = new BasketCounter_t[buffer_size];
    fBufferMask = buffer_size - 1;
    fBsize = basket_size;
    fBmask = ~(fBsize - 1);
    for (size_t i = 0; i < buffer_size; ++i) {
      fBuffer[i] = nullptr;
      fCounters[i].SetBsize(fBsize);
    }
  }

  //____________________________________________________________________________
  /** @brief Add an element to the basketizer */
  GEANT_INLINE
  bool AddElement(T *const data, std::vector<T *> &basket) {
    // Book atomically a slot for copying the element
    size_t ibook = fIbook.fetch_add(1);
    // Compute position in the buffer
    size_t ibook_buf = ibook & fBufferMask;
    // Get basket address for the booked index. Assume fixed size baskets.
    size_t ibasket = ibook_buf & fBmask;
    size_t nbooked = fCounters[ibasket].BookSlot();
    if (ibasket == ibook_buf)
      fNbaskets++;
    size_t nfilled = fCounters[ibasket].FillSlot(nbooked, &fBuffer[ibook_buf], data);
    fNstored++;
    assert (nfilled <= fBsize);
    if (nfilled == fBsize) {
      // Deploy basket (copy data)
      basket.insert(basket.end(), &fBuffer[ibasket], &fBuffer[ibasket + fBsize]);
      fNstored -= fBsize;
      fNbaskets--;
      fCounters[ibasket].ReleaseBasket();
      return true;
    }
    return false;
  }

  //____________________________________________________________________________
  /** @brief Garbage collect data */
  GEANT_INLINE
  bool GarbageCollect(std::vector<T *> &basket) {
    // Garbage collect all elements present in the container on the current basket
    Lock();
    // Get current index, then compute a safe location forward
    size_t ibook = fIbook.load(std::memory_order_acquire);
    size_t inext = (ibook + 1 + fBsize) & fBufferMask;
    size_t inext_start = inext & fBmask;
    while (!fIbook.compare_exchange_weak(ibook, inext_start, std::memory_order_relaxed)) {
      ibook = fIbook.load(std::memory_order_acquire);
      inext = (ibook + 1 + fBsize) & fBufferMask;
      inext_start = inext & fBmask;
    }
    // Now we can cleanup the basket at ibook
    size_t buf_pos = ibook & fBufferMask;
    size_t buf_start = buf_pos & fBmask;
    size_t nelem = buf_pos - buf_start;
    if (!nelem) {
      Unlock();
      return false;
    }
    while (fCounters[buf_start].Nfilled() < nelem)
      ;
    fNstored -= nelem;
    // Now fill the objects from the pending basket
    basket.insert(basket.end(), &fBuffer[buf_start], &fBuffer[buf_start + nelem]);
    fNbaskets--;
    fCounters[buf_start].Clear();
    Unlock();
    return true;
  }

  //____________________________________________________________________________
  /** @brief Change dynamically basket size */
  GEANT_INLINE
  bool SetBasketSize(unsigned int bsize, std::vector<T *> & /*basket*/) {
    int shift = 0;
    size_t basket_size;
    while (bsize >> (++shift + 1))
      ;
    basket_size = 1 << shift;
    if (basket_size == fBsize)
      return false;
    // Do not allow any other garbage collections or changing of basket size
    // while doing this.
    Lock();
    // The algorithm here will need to deal with usage of the previous fBsize
    // concurrently by AddElement
    Unlock();
    return false;
  }

  //____________________________________________________________________________
  /** @brief GetNumber of pending baskets */
  GEANT_INLINE
  short int GetNpending() const { return fNbaskets.load(); }

  //____________________________________________________________________________
  /** @brief Lock the container for GC */
  GEANT_INLINE
  void Lock() {
    while (fLock.test_and_set(std::memory_order_acquire))
      ;
  }

  //____________________________________________________________________________
  /** @brief Unlock the container for GC */
  GEANT_INLINE
  void Unlock() { fLock.clear(std::memory_order_release); }

private:
  static const size_t cacheline_size = 64;
  typedef char cacheline_pad_t[cacheline_size];
  size_t fBsize;          /** Size of the produced baskets */
  size_t fBmask;          /** Basket mask */
  std::atomic_flag fLock; /** Lock for garbage collection and bsize changes */
  // Make sure the following data members are not sitting in the same cache line
  cacheline_pad_t pad0_;
  T **fBuffer; /** Circular buffer for elements to be basketized*/
  cacheline_pad_t pad1_;
  BasketCounter_t *fCounters; /** Basket counters */
  cacheline_pad_t pad2_;
  size_t fBufferMask; /** Buffer mask used for fast index calculation in buffer */
  cacheline_pad_t pad3_;
  std::atomic<size_t> fIbook; /** Current booked index */
  cacheline_pad_t pad4_;
  std::atomic<size_t> fNstored; /** Number of objects currently stored */
  cacheline_pad_t pad5_;
  std::atomic<short> fNbaskets; /** Number of pending baskets */
  cacheline_pad_t pad6_;
};
} // Geant

#endif
