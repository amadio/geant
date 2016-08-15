//===--- GeantVTaskMgr.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantVTaskMgr.h
 * @brief Task management interface.
 * @author andrei.gheata@cern.ch
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_VTASKMGR
#define GEANT_VTASKMGR

class WorkloadManager;

/** @brief GeantVTaskMgr class */
class GeantVTaskMgr {
public:
  
  /** @brief GeantVTaskMgr constructor */	
  GeantVTaskMgr() {}

  /** @brief GeantVTaskMgr destructor */
  virtual ~GeantVTaskMgr() {}

  /** @brief Function of initialization */
  virtual bool Initialize(int nthreads) = 0;

};
#endif
