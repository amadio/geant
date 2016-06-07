//===--- NumaBlock.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaBlock.h
 * @brief NUMA-aware data blocks
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_BLOCK
#define GEANT_NUMA_BLOCK

#include <atomic>
#include <cassert>
#include <iostream>
#include <type_traits>
#include "Geant/Config.h"

namespace Geant {
/**
 * @brief Class NumaBlock
 * @detailed A NUMA block is a concurrent templated factory holding contiguously a
 *           number of objects. It keeps a counters for numbers of used objects.
 */
template <typename T> class NumaBlock {

  using size_t = std::size_t size_t;
  using atomic_size_t = std::atomic<std::size_t>;
  static size_t const cacheline_size = 64;
  typedef char cacheline_pad_t[cacheline_size];

private:
  atomic_size_t fCurrent;  // Current free track index
  atomic_size_t fUsed;     // Number of tracks used

  cacheline_pad_t pad0_;   //! Padding to protect the other data from the hot cache line above

  size_t        fSize;     // Number of tracks stored by the block
  T             fArray[1]; //! Array of elements
  

private:
  /** @brief Constructor */
  NumaBlock(size_t size) : fCurrent(0), fUsed(0), fSize(size), fArray(0)
  {
    // NUMA block constructor. If the system is NUMA-aware, the block will be alocated
    // on the memory associated with the given NUMA node.
    static_assert(!std::is_polymorphic<T>::value, "Cannot use polymorphic types as GeantBlock");
    static_assert(std::is_default_constructible<T>::value, "Type used in GeantBlock must have default ctor.");
    static_assert(std::is_copy_constructible<T>::value, "Type used in GeantBlock must be copy constructible");

    fArray = reinterpret_cast<T*>(numa_aligned_malloc(fSize * sizeof(T), numa_node, 64));
    for (size_t i=0; i<size; ++i) new (&fArray[i]) T();
  }
  
  NumaBlock(const NumaBlock&);
  NumaBlock& operator=(const NumaBlock&);

public:
  static NumaBlock *MakeInstance(size_t nvalues, int numa_node)
  {
    // Make an instance. To be released using ReleaseInstance. 
    size_t needed = SizeOf(nvalues);
    void *ptr = numa_aligned_malloc(needed, numa_node, 64);
    NumaBlock *block = new (ptr) NumaBlock(nvalues);
    return ( block );    
  }
  
  static constexpr size_t SizeOf(size_t nvalues)
  { return ( sizeof(NumaBlock<T>) + nvalues*sizeof(T) ); }
  
public:
  GEANT_INLINE T *GetValues() { return &fArray[0]; }
  GEANT_INLINE const T *GetValues() const { return &fArray[0]; }

  GEANT_INLINE T &operator[](size_t index) { return GetValues()[index]; };
  GEANT_INLINE const T &operator[](size_t index) const { return GetValues()[index]; };
  

  /** @brief Destructor */
  ~TrackBlock() { numa_aligned_free(fArray); }
  
  /** @brief Get an object pointer from the container */
  GEANT_INLINE T *GetObject() {
    size_t current = fCurrent.fetch_add(1);
    if (current >= fSize) return nullptr;
    fUsed++;
    return ( &fArray[current] );
  }
  
  /** @brief Release an object to the container */
  GEANT_INLINE size_t ReleaseObject() { return ( fUsed.fetch_sub(1) - 1 ); }
  
  /** @brief Check if the block is still in use */
  GEANT_INLINE bool InUse() const { return (fUsed.load() == 0); }
  /** @brief Check if track belongs to container */
  GEANT_INLINE bool OwnsObject(T *obj) const {
    return ( (size_t)(track - fArray) < fSize*sizeof(T) );
  }
};

} // Geant
