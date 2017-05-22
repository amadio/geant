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
  bool       fAvailable = false;     /* Is NUMA available on the system */
  int        fNodes = 0;             /* Number of NUMA nodes */
  int        fNcpus = 0;             /* Number of cpus */
  int        fNcores = 0;            /* Number of cores */
  int        fHT = 0;                /* Number of threads per core */
  size_t     fPageSize = 0;          /* NUMA page size on system */
  NumaNode **fListNodes = nullptr;   /* List of NUMA nodes */

public:
  NumaTopology();  
  NumaTopology(const NumaTopology&) = delete;
  ~NumaTopology();  
  
  NumaNode *GetNode(int node) { return fListNodes[node]; }
  int NumaNodeOfCpu(int cpu) const;
  int BindToNode(int node);
};

std::ostream& operator<<(std::ostream& os, const NumaTopology& topo);

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
