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

/* NUMA aligned memory allocator */
void *numa_aligned_malloc(std::size_t bytes, int node, std::size_t alignment);
void  numa_aligned_free(void *p);

/* NUMA memory address inspector */
int   numa_node_addr(void *ptr);

/* Pin a thread to a core */
void pin_to_core(size_t core);

} // Geant

#endif
