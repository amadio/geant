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

#include "Geant/Config.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantPropagator;

/** @brief GeantVTaskMgr class */
class GeantVTaskMgr {
public:
  
  /** @brief GeantVTaskMgr constructor */	
  GeantVTaskMgr() {}

  /** @brief GeantVTaskMgr destructor */
  virtual ~GeantVTaskMgr() {}

  /** @brief Function of initialization */
  virtual bool Initialize(int nthreads, GeantPropagator *prop) = 0;

  /** @brief Function for final actions */
  virtual void Finalize() = 0;
};

} // GEANT_IMPL_NAMESPACE
} // Geant
#endif
