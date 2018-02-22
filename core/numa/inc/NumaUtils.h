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

#ifdef GEANT_USE_NUMA
#include <mutex>
#include <hwloc.h>
#endif

#include "Geant/Config.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

namespace NumaUtils {
  VECCORE_ATT_HOST_DEVICE
  void *NumaAlignedMalloc(std::size_t bytes, int node, std::size_t alignment);

  VECCORE_ATT_HOST_DEVICE
  void  NumaAlignedFree(void *p);

  /* @brief NUMA memory address inspector */
  VECCORE_ATT_HOST_DEVICE
  int   NumaNodeAddr(void *ptr);
  
  /* @brief Pin a thread to a core */
  VECCORE_ATT_HOST_DEVICE
  int GetCpuBinding();
  
  VECCORE_ATT_HOST_DEVICE
  bool NumaAvailable();

#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
  hwloc_topology_t const &Topology();
#endif
}

#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)

struct NumaUtilsStruct {
  static NumaUtilsStruct *fgInstance; /** Singleton instance */
  bool fAvailable = false;
  std::mutex fLock;
  hwloc_topology_t fTopology; /* NUMA topology context */

  /** @brief Constructor **/
  VECCORE_ATT_HOST_DEVICE
  NumaUtilsStruct();

  /** @brief Destructor **/
  VECCORE_ATT_HOST_DEVICE
  ~NumaUtilsStruct();
  
  /** @brief Function that creates NumaUtils instance **/
  VECCORE_ATT_HOST_DEVICE
  static NumaUtilsStruct* Instance();

  
  /* @brief NUMA aligned memory allocator */
  VECCORE_ATT_HOST_DEVICE
  void *NumaAlignedMalloc(std::size_t bytes, int node, std::size_t alignment);

  VECCORE_ATT_HOST_DEVICE
  void  NumaAlignedFree(void *p);

  /* @brief NUMA memory address inspector */
  VECCORE_ATT_HOST_DEVICE
  int   NumaNodeAddr(void *ptr);
  
  /* @brief Pin a thread to a core */
  VECCORE_ATT_HOST_DEVICE
  int GetCpuBinding() const;
  
};

#endif

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
