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

namespace Geant {

/** Class describung NUMA node properties */
//______________________________________________________________________________
class NumaNode {
public:
  int        fId;            /* NUMA node id */
  int        fNcpus;         /* Number of CPU's allowed for the node */
  int        fTotcpus;       /* Total number of CPU's */
  int        fNphysical;     /* Number of physical cores */
  int        fNthreads;      /* Number of assigned threads for the node */ 
  long       fMemTotal;      /* Total memory */
  long       fMemFree;       /* Free memoru */
  int       *fCpus;          /* List of CPU's for the node */
  std::mutex fMutex;         /* Mutex for the node */
  
public:
  NumaNode(int id, int maxcpus);  
  ~NumaNode() { delete [] fCpus; }
  
  int PinThread();
};

std::ostream& operator<<(std::ostream& os, const NumaNode& node);

} // Geant

#endif
