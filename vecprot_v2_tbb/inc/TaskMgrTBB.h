//===--- TaskMgrTBB.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file TaskMgrTBB.h
 * @brief TBB task manager.
 * @author andrei.gheata@cern.ch
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TASKMGRTBB
#define GEANT_TASKMGRTBB

#include "GeantVTaskMgr.h"

class WorkloadManager;

/** @brief TaskMgrTBB class */
class TaskMgrTBB : public GeantVTaskMgr {
public:
  
  /** @brief TaskMgrTBB constructor */	
  TaskMgrTBB() {}

  /** @brief TaskMgrTBB destructor */
  virtual ~TaskMgrTBB() {}

  /** @brief Function of initialization */
  virtual bool Initialize(int nthreads);

};
#endif
