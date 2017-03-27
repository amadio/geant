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
#include <cassert>
#include <iostream>
#include "Geant/Typedefs.h"
#include "Geant/Config.h"


/**
 * @brief Class Basketizer
 *
 * @param buffer_size Buffer size for queue
 * @tparam T Type of objects
 */
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

template <typename O, typename I>
VECCORE_ATT_HOST_DEVICE
O fetch_add(O &val, I add) {
  auto old = val;
  val += add;
  return old;
}

template <typename O, typename I>
VECCORE_ATT_HOST_DEVICE
O fetch_add_ret_newval(O &val, I add) {
  auto old = val;
  val += add;
  return old;
}

template <typename O, typename I>
O fetch_add(std::atomic<O> &val, I add) {
  return (val.fetch_add(add));
}

template <typename O, typename I>
O fetch_add_ret_newval(std::atomic<O> &val, I add) {
  return (val.fetch_add(1) + add);
}

template <typename T> class BasketCounter {
  using size_t = std::size_t;

#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  // On cuda there is one propagator per thread.  So (for now), no need
  // for atomics.
  template <typename N> using atomic_t = N;
#else
  template <typename N>
  using atomic_t = std::atomic<N>;
#endif

public:
  size_t fBsize;                // Basket size
  atomic_t<size_t> fIbook;   // Booking base index
  atomic_t<size_t> fNbook0;  // Booking base counter
  atomic_t<size_t> fNbooktot; // Counter for booked slots
  atomic_t<size_t> fNfilled; // Counter for filled slots

  VECCORE_ATT_HOST_DEVICE
  BasketCounter() : fBsize(0), fNbooktot(0), fNfilled(0) { }

  VECCORE_ATT_HOST_DEVICE
  BasketCounter(short bsize) : fBsize(bsize), fIbook(0), fNbook0(0), fNbooktot(0), fNfilled(0) { }

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetBsize(size_t bsize) { fBsize = bsize; }

  //____________________________________________________________________________
  GEANT_FORCE_INLINE
  size_t Bsize() { return (fBsize); }

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t BookSlot(size_t ibook) {
    while (ibook - Ibook() >= fBsize)
      ;
    return fetch_add_ret_newval(fNbooktot,1); // fNbooktot.fetch_add(1) + 1);
  }

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool BookSlots(size_t ibook, size_t expected, size_t nslots) {
     if (ibook - Ibook() >= fBsize)
       return false;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
     while (!fNbooktot.compare_exchange_weak(expected, expected+nslots, std::memory_order_relaxed))
       ;
#else
    fNbooktot = expected+nslots;
#endif
     return true;
  }

  //____________________________________________________________________________
  void Print() {
    std::cout << "== ibook = " << Ibook() << " nbook0 = " << Nbook0()
              << " nbooktot = " << Nbooktot() << " nfilled = " << Nfilled() << std::endl;
  }

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t FillSlot(size_t nbooktot, T *address, T const data) {
    // Block thread until booking range matches the current basket
    // This can create contention if the buffer is too small or the
    // workload is too small.
    while (nbooktot - Nbook0() > fBsize)
      ;
    // Copy data to the booked slot
    *address = data;
    return fetch_add_ret_newval(fNfilled,1); // (fNfilled.fetch_add(1) + 1);
  }

  //____________________________________________________________________________
  GEANT_FORCE_INLINE
  size_t FillDummy(size_t nslots) {
    return fetch_add_ret_newval(fNfilled,nslots); // (fNfilled.fetch_add(nslots) + nslots);
  }  

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void ReleaseBasket(size_t ibook) {
    fNfilled = 0;
    fNbook0 += fBsize;
    fIbook = ibook;
  }

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetIbook(size_t ibook) { fIbook = ibook; }

  //____________________________________________________________________________
  GEANT_FORCE_INLINE
  void Clear(size_t nclear) {
    fNfilled -= nclear;
    fNbook0 += nclear;
  }

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t Ibook() const { return fIbook; }

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t Nbooktot() const { return fNbooktot; }

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t Nbook0() const { return fNbook0; }

  //____________________________________________________________________________
  // Note: to calculate the result, this does an arithmetic operation on two
  // atomic but do not guarantee that they are fetch at the exact same time.
  // I.e. the result is approximatif.
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t Nbooked() const { return fNbooktot - fNbook0; }

  //____________________________________________________________________________
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t Nfilled() const { return fNfilled; }
};

template <typename T> class Basketizer {
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  // On cuda there is one propagator per thread.  So (for now), no need
  // for atomics.
  template <typename N> using atomic_t = N;
#else
  template <typename N>
  using atomic_t = std::atomic<N>;
#endif
public:
  using size_t = std::size_t;
  using BasketCounter_t = BasketCounter<T *>;

  /**
   * @brief Basketizer dummy constructor
   */
  Basketizer()
      : fBsize(0), fBmask(0), fBuffer(0), fCounters(0), fBufferMask(0), fIbook(0), fNstored(0), fNbaskets(0) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    fLock.clear();
#endif
  }


  /**
   * @brief Basketizer default constructor
   *
   * @param buffer_size Circular buffer size
   * @param basket_size Size of produced baskets
   */
  VECCORE_ATT_HOST_DEVICE
  Basketizer(size_t buffer_size, unsigned int basket_size, void *addr=0)
      : fBsize(0), fBmask(0), fBuffer(0), fCounters(0), fBufferMask(0), fIbook(0), fNstored(0), fNbaskets(0) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    fLock.clear();
#endif
    Init(buffer_size, basket_size, addr);
  }

  /** @brief Basketizer destructor */
  VECCORE_ATT_HOST_DEVICE
  ~Basketizer() {
    delete[] fBuffer;
    delete[] fCounters;
  }

  //____________________________________________________________________________
  /** @brief Initialize basketizer */
  VECCORE_ATT_HOST_DEVICE
  void Init(size_t buffer_size, unsigned int basket_size, void *addr = 0) {
    // Make sure the requested size is a power of 2
    assert((basket_size >= 2) && ((basket_size & (basket_size - 1)) == 0));
    fBsize = basket_size;
    if (addr) {
      fBuffer = reinterpret_cast<T**>((char*)addr + sizeof(Basketizer<T>));
      fCounters = reinterpret_cast<BasketCounter_t*>((char*)addr + sizeof(Basketizer<T>) + buffer_size * sizeof(T*));
    } else {
      fBuffer = new T *[buffer_size];
      fCounters = new BasketCounter_t[buffer_size];
    }
    fBufferMask = buffer_size - 1;
    fBsize = basket_size;
    fBmask = ~(fBsize - 1);
    for (size_t i = 0; i < buffer_size; ++i) {
      fBuffer[i] = nullptr;
      if (i%basket_size == 0) {
        fCounters[i].SetBsize(fBsize);
        fCounters[i].SetIbook(i);
      }  
    }
  }

  //____________________________________________________________________________
  /** @brief Add an element to the basketizer */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool AddElement(T *const data, vector_t<T *> &basket) {
    // Book atomically a slot for copying the element
    size_t ibook = fetch_add(fIbook,1);
    // Compute position in the buffer
    size_t ibook_buf = ibook & fBufferMask;
    // Get basket address for the booked index. Assume fixed size baskets.
    size_t ibasket = ibook_buf & fBmask;
    size_t nbooktot = fCounters[ibasket].BookSlot(ibook);
    if (ibasket == ibook_buf)
      fNbaskets++;
    size_t nfilled = fCounters[ibasket].FillSlot(nbooktot, &fBuffer[ibook_buf], data);
    fNstored++;
    if (nfilled == fBsize) {
      // Deploy basket (copy data)
#ifndef VECCORE_CUDA
      basket.insert(basket.end(), &fBuffer[ibasket], &fBuffer[ibasket + fBsize]);
#else
      auto end = ibasket + fBsize;
      for(auto cursor = ibasket; cursor < end; ++cursor) {
        basket.push_back(fBuffer[cursor]);
      }
#endif
      fNstored -= fBsize;
      fNbaskets--;
      fCounters[ibasket].ReleaseBasket(fBufferMask + 1 + (ibook & fBmask));
      return true;
    }
    return false;
  }

  //____________________________________________________________________________
  /** @brief Check baskets for remaining tracks */
  GEANT_FORCE_INLINE
  void CheckBaskets() {
    size_t nbaskets = (fBufferMask+1)/fBsize;
    Lock();
    std::cout << "fNstored = " << fNstored << "  fNbaskets = " << fNbaskets << std::endl;
    size_t nstored = 0;
    size_t nbooked = 0;
    size_t nbooktot = 0;
    for (size_t ib=0; ib<nbaskets; ++ib) {
      size_t ibasket = ib*fBsize;
      nstored += fCounters[ibasket].Nfilled();
      nbooktot += fCounters[ibasket].Nbooktot();
      nbooked += fCounters[ibasket].Nbooked();
    }
    std::cout << "ibook = " << fIbook.load() << "  sumbooked = " << nbooktot << "  sumstored = " << nstored << "  stillbooked = " << nbooked << std::endl;  
    Unlock();
  }  

  //____________________________________________________________________________
  /** @brief Garbage collect data */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool Flush(vector_t<T *> &basket) {
    // Flush all elements present in the container on the current basket
    // If this is not the last basket, drop garbage collection
    if (fNbaskets > 1)
      return false;
    Lock();
    // Get current index, then compute a safe location forward
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    size_t ibook = fIbook.load(std::memory_order_acquire);
#else
    size_t ibook = fIbook;
#endif
    size_t inext = (ibook + fBsize) & fBmask;
    // Now we can cleanup the basket at ibook.
    size_t ibook_buf = ibook & fBufferMask;
    size_t ibasket = ibook_buf & fBmask;
    size_t nelem = ibook_buf - ibasket;
    size_t nbooktot = fCounters[ibasket].Nbooktot();
    size_t nbooked = fCounters[ibasket].Nbooked();
    size_t nfilled = fCounters[ibasket].Nfilled();
    size_t ibook0 = fCounters[ibasket].Ibook();
    if ((!nelem) || (nbooked!=nelem) || (nfilled != nelem) || (ibook-ibook0>=fBsize)) {
      Unlock();
      return false;
    }
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    if (!fIbook.compare_exchange_strong(ibook, inext, std::memory_order_relaxed)) {
      Unlock();
      return false;    
    }
#else
    fIbook = inext;
#endif
    // We managed to push ibook to a safe location, so new bookers cannot book this
    // basket till it gets released. Now book the remaining slots in the basket
    while (!fCounters[ibasket].BookSlots(ibook, nbooktot, fBsize-nelem))
      ;
    // Now fill the objects from the pending basket
#ifndef VECCORE_CUDA
    basket.insert(basket.end(), &fBuffer[ibasket], &fBuffer[ibasket + nelem]);
#else
    auto end = ibasket + fBsize;
      for(auto cursor = ibasket; cursor < end; ++cursor) {
        basket.push_back(fBuffer[cursor]);
      }
#endif

    fNstored -= nelem;
    fNbaskets--;
    fCounters[ibasket].ReleaseBasket(fBufferMask + 1 + (ibook & fBmask));
    
    Unlock();
    return true;
  }

  //____________________________________________________________________________
  /** @brief Change dynamically basket size */
  GEANT_FORCE_INLINE
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
  GEANT_FORCE_INLINE
  short int GetNpending() const { return fNbaskets.load(); }

  //____________________________________________________________________________
  /** @brief GetNumber of stored objects */
  GEANT_FORCE_INLINE
  size_t GetNstored() const { return fNstored.load(); }

  //____________________________________________________________________________
  /** @brief Lock the container for GC */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void Lock() {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    while (fLock.test_and_set(std::memory_order_acquire))
      ;
#endif
  }

  //____________________________________________________________________________
  /** @brief Unlock the container for GC */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void Unlock() {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    fLock.clear(std::memory_order_release);
#endif
  }

  //____________________________________________________________________________
  /** @brief Get size of a basketizer instance depending on the buffer size */
  static
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t SizeofInstance(size_t buffer_size) {
    return (sizeof(Basketizer<T>) + buffer_size * (sizeof(T*) + sizeof(BasketCounter_t)));
  }
  
  //____________________________________________________________________________
  /** @brief Make basketizer instance at a given address */
  static
  VECCORE_ATT_HOST_DEVICE
  Basketizer<T> *MakeInstanceAt(void *addr, size_t buffer_size, size_t basket_size) {
    return new (addr) Basketizer<T>(buffer_size, basket_size);
  }

private:
  static const size_t cacheline_size = 64;
  typedef char cacheline_pad_t[cacheline_size];
  size_t fBsize;          /** Size of the produced baskets */
  size_t fBmask;          /** Basket mask */
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  std::atomic_flag fLock; /** Lock for garbage collection and bsize changes */
#endif
    // Make sure the following data members are not sitting in the same cache line
  cacheline_pad_t pad0_;
  T **fBuffer; /** Circular buffer for elements to be basketized*/
  cacheline_pad_t pad1_;
  BasketCounter_t *fCounters; /** Basket counters */
  cacheline_pad_t pad2_;
  size_t fBufferMask; /** Buffer mask used for fast index calculation in buffer */
  cacheline_pad_t pad3_;
  atomic_t<size_t> fIbook; /** Current booked index */
  cacheline_pad_t pad4_;
  atomic_t<size_t> fNstored; /** Number of objects currently stored */
  cacheline_pad_t pad5_;
  atomic_t<short> fNbaskets; /** Number of pending baskets */
  cacheline_pad_t pad6_;
};
} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
