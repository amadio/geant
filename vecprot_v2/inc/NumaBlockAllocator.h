//===--- NumaBlockAllocator.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaBlockAllocator.h
 * @brief Memory manager for NUMA-aware data blocks
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_BLOCK_ALLOCATOR
#define GEANT_NUMA_BLOCK_ALLOCATOR

#include <atomic>
#include <vector>
#include <cassert>
#include <iostream>
#include <type_traits>
#include "Geant/Config.h"

/**
 * @brief Class NumaBlockAllocator
 * @detailed The class is managing allocation of blocks of POD data for a given 
 *           NUMA node. The allocation for single blocks is using numa_aligned_malloc.
 */
namespace Geant {
/**
 * @brief Class NumaBlock
 * @detailed A NUMA block is a concurrent templated factory holding contiguously a
 *           number of objects. It keeps a counters for numbers of used objects.
 */
template <typename T> class NumaBlock {
  using size_t = std::size_t size_t;
  using atomic_size_t = std::atomic<std::size_t>;
private:
  size_t fSize;           // Number of tracks stored by the block
  atomic_size_t fCurrent; // Current free track index
  atomic_size_t fUsed;    // Number of tracks used
  T*            fArray;   // Array of elements
  
public:
  /** @brief Constructor */
  NumaBlock(size_t size, int numa_node=0) : fSize(size), fCurrent(0), fUsed(0), fTracks(0)
  {
    // NUMA block constructor. If the system is NUMA-aware, the block will be alocated
    // on the memory associated with the given NUMA node.
    static_assert(!std::is_polymorphic<T>::value, "Cannot use polymorphic types as GeantBlock");
    static_assert(std::is_default_constructible<T>::value, "Type used in GeantBlock must have default ctor.");
    static_assert(std::is_copy_constructible<T>::value, "Type used in GeantBlock must be copy constructible");

    fArray = reinterpret_cast<T*>(numa_aligned_malloc(fSize * sizeof(T), numa_node, 64));
    for (size_t i=0; i<size; ++i) new (&fArray[i]) T();
  }


  /** @brief Destructor */
  ~TrackBlock() { numa_aligned_free(fArray); }
  
  /** @brief Get an object pointer from the container */
  GEANT_INLINE
  T *GetObject() {
    size_t current = fCurrent.fetch_add(1);
    if (current >= fSize) return nullptr;
    fUsed++;
    return ( &fArray[current] );
  }
  
  /** @brief Release an object to the container */
  GEANT_INLINE
  size_t ReleaseObject() { return ( fUsed.fetch_sub(1) - 1 ); }
  
  /** @brief Check if track belongs to container */
  GEANT_INLINE
  bool OwnsObject(T *obj) const {
    return ( (size_t)(track - fArray) < fSize*sizeof(T) );
  }
    

};
} // Geant
