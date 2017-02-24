//===--- DiscreteProcStage.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file DiscreteProcStage.h
 * @brief The stage dealing with discrete processes.
 * @details Applies discrete processes to all staged tracks
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_discrete_PROC_STAGE
#define GEANT_discrete_PROC_STAGE

#include "SimulationStage.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class DiscreteProcStage : public SimulationStage {

private:
  DiscreteProcStage(const DiscreteProcStage &) = delete;
  DiscreteProcStage &operator=(const DiscreteProcStage &) = delete;

protected:

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

public:

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(GeantTrack *track);

public:
  /** @brief Dummy DiscreteProcStage constructor */
  VECCORE_ATT_HOST_DEVICE
  DiscreteProcStage() {}

  /** @brief Standard DiscreteProcStage constructor */
  VECCORE_ATT_HOST_DEVICE
  DiscreteProcStage(GeantPropagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  ~DiscreteProcStage();

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
