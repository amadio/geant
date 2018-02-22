//===--- NumaPolicy.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaPolicy.h
 * @brief Class for implementing thread allocation NUMA policies
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_POLICY
#define GEANT_NUMA_POLICY

#include "Geant/Config.h"

#ifndef GEANT_NUMA_TOPOLOGY
#include "NumaTopology.h"
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

/* NUMA policy class  */
//______________________________________________________________________________
class NumaPolicy {
public:
  enum EPolicyType {
    kSysDefault = 0,       /** System default scheduling policy */
    kCompact    = 1 << 0,  /** Compact allocation until filling every NUMA node */
    kScatter    = 1 << 1,  /** Scatter threads evenly across the system */
    kHTcompact  = 1 << 2   /** Use HT in compact mode before pinning to next NUMA node */
  };

  int          fNthreads;        /* Number of threads to be pinned */
  NumaTopology fTopo;            /* Numa topology of the machine */
  EPolicyType  fPolicy;          /* NUMA policy */
  
  NumaPolicy(EPolicyType policy) : fNthreads(0), fTopo(), fPolicy(policy) {}
  
  NumaTopology *GetTopology() { return &fTopo; }
  void SetPolicy(EPolicyType policy) { fPolicy = policy; }

  int AllocateNextThread(int node = -1);
  int MembindNode(int node);
  int NextNumaNode();
  int GetNnumaNodes() { return fTopo.fNodes; }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
