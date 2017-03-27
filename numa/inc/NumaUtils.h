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

#ifndef GEANT_NUMA_TOPOLOGY
#include "NumaTopology.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/* NUMA aligned memory allocator */
VECCORE_ATT_HOST_DEVICE
void *NumaAlignedMalloc(std::size_t bytes, int node, std::size_t alignment);

VECCORE_ATT_HOST_DEVICE
void  NumaAlignedFree(void *p);

/* NUMA memory address inspector */
int   NumaNodeAddr(void *ptr);

/* Pin a thread to a core */
void PinToCore(size_t core);

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
