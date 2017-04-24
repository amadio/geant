//===--- NumaTopology.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaTopology.h
 * @brief Class describing the NUMA topology of the system, allowing to set NUMA policies
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_TOPOLOGY
#define GEANT_NUMA_TOPOLOGY

#include <iostream>
#include "Geant/Config.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class NumaNode;

/* NUMA topology class storing NUMA-specific data detected for the system */
//______________________________________________________________________________
class NumaTopology {
public:
  bool       fAvailable;     /* Is NUMA available on the system */
  int        fNodes;         /* Number of allowed NUMA nodes */
  int        fNcpus;         /* Number of allowed CPU's */
  int        fNphysical;     /* Number of physical cores */
  int        fHT;            /* Number of threads per core */
  size_t     fPageSize;      /* NUMA page size on system */
  NumaNode **fListNodes;     /* List of nodes */
  int       *fNthreads;      /* Number of threads currently allocated per NUMA node */

public:
  NumaTopology();  
  NumaTopology(const NumaTopology&) = delete;
  ~NumaTopology();  
  
  int FindCores();
  NumaNode *GetNode(int node) { return fListNodes[node]; }
  
  int PinToNode(int node);
};

std::ostream& operator<<(std::ostream& os, const NumaTopology& topo);

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
