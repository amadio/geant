//===--- NumaNode.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaNode.h
 * @brief Class for implementing a NUMA node
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_NODE
#define GEANT_NUMA_NODE

#include <mutex>
#include <iostream>
#include "Geant/Config.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class NumaCore;

/** Class describung NUMA node properties */
//______________________________________________________________________________
class NumaNode {
public:
  int        fId = -1;         /* NUMA node id */
  int        fNcores = 0;      /* Number of physical cores */
  int        fNcpus = 0;       /* Number of CPU's on this node */
  int        fNthreads = 0;    /* Number of assigned threads for the node */ 
  long       fMemTotal = 0;    /* Total memory */
  long       fMemFree = 0;     /* Free memory */
  NumaCore **fCores = nullptr; /* List of CPU's for the node */
  std::mutex fMutex;           /* Mutex for the node */
  
public:
  NumaNode(int id);
  ~NumaNode();
  NumaNode(const NumaNode&) = delete;
  NumaNode &operator=(const NumaNode&) = delete;
  
  bool HasCpu(int cpu) const;  
  int BindThread();
};

std::ostream& operator<<(std::ostream& os, const NumaNode& node);

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
