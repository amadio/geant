//===--- NumaUtils.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaUtils.h
 * @brief Utilities for NUMA management
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_UTILS
#define GEANT_NUMA_UTILS

#ifdef USE_NUMA
#include <mutex>
#include <hwloc.h>
#endif

#include "Geant/Config.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

struct NumaUtils {
  static NumaUtils *fgInstance; /** Singleton instance */
  bool fAvailable = false;
#ifdef USE_NUMA
  std::mutex fLock;
  hwloc_topology_t fTopology; /* NUMA topology context */
#endif

  /** @brief Constructor **/
  VECCORE_ATT_HOST_DEVICE
  NumaUtils();

  /** @brief Destructor **/
  VECCORE_ATT_HOST_DEVICE
  ~NumaUtils();
  
  /** @brief Function that creates NumaUtils instance **/
  static NumaUtils* Instance();

  
  /* @brief NUMA aligned memory allocator */
  VECCORE_ATT_HOST_DEVICE
  void *NumaAlignedMalloc(std::size_t bytes, int node, std::size_t alignment);

  VECCORE_ATT_HOST_DEVICE
  void  NumaAlignedFree(void *p);

  /* @brief NUMA memory address inspector */
  int   NumaNodeAddr(void *ptr);
  
  /* @brief Pin a thread to a core */
  int GetCpuBinding() const;
  
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
