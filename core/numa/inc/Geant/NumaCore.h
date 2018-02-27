//===--- NumaCore.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaCore.h
 * @brief Class representing a processor core
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_CORE
#define GEANT_NUMA_CORE

#include <iostream>
#include "Geant/NumaUtils.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
class NumaCore {
public:
  int         fId = -1;           /* NUMA core id */
  int         fNode = -1;         /* NUMA node where it belongs */
  int         fNcpus = 0;         /* Number of logical CPU's for the core */
  int         fNthreads = 0;      /* Number of assigned threads for the core */
  int        *fCpus = nullptr;    /* List of CPU's for the node */
#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
  std::mutex  fMutex;             /* Mutex for the node */
  hwloc_obj_t fObjCore = nullptr; /* Object in the topology representing this core */
#endif
public:
  NumaCore(int id, int node);
  ~NumaCore() { delete [] fCpus; }
  
  bool HasCpu(int cpu) const;  
  int BindThread();
};

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& os, const NumaCore& core);

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
