//===--- NumaAllocator.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaAllocator.h
 * @brief A NUMA aware allocator
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_ALLOCATOR
#define GEANT_NUMA_ALLOCATOR

#include <iostream>

#ifndef GEANT_NUMA_UTILS
#include "Geant/NumaUtils.h"
#endif

/**
  NUMA-aware allocator following std::alocator requirements.
  Use with any type as:
    std::vector<T, NumaAllocator<T>> basket(init_size, init_value, NumaAllocator<T>(numafNode));
*/

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

template <typename T>
class NumaAllocator : public std::allocator<T> {
private:
  int fNode; /** Numa node */

public:
  typedef size_t size_type;
  typedef T *pointer;
  typedef const T *const_pointer;

  template <typename U>
  struct rebind {
    typedef NumaAllocator<U> other;
  };

  pointer allocate(size_type n)
  {
    // std::cout << "Alloc " << n*sizeof(T) << " bytes on numa node " << fNode << std::endl;
    return static_cast<pointer>(NumaUtils::NumaAlignedMalloc(n * sizeof(T), fNode, 64));
    // return std::allocator<T>::allocate(n, hint);
  }

  void deallocate(pointer p, size_type /*n*/)
  {
    // std::cout <<  "Dealloc " <<  n*sizeof(T) << " bytes.\n";
    NumaUtils::NumaAlignedFree(p);
    // return std::allocator<T>::deallocate(p, n);
  }

  NumaAllocator() throw() : std::allocator<T>(), fNode(0) {}
  NumaAllocator(int node) throw() : std::allocator<T>(), fNode(node) {}
  NumaAllocator(const NumaAllocator &a) throw() : std::allocator<T>(a), fNode(a.fNode) {}

  template <class V>
  NumaAllocator(const NumaAllocator<V> &a) throw() : std::allocator<T>(a), fNode(a.fNode)
  {
  }

  ~NumaAllocator() throw() {}
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
