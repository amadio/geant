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
    
  /** @brief Retrieve the handler for vectorized processes, or the default handler
   *  === TO BE IMPLEMENTED FOR ALREADY VECTORIZED PROCESSES ===
   */
  VECCORE_ATT_HOST_DEVICE
  Handler *GetHandler(int /*process*/) { return fHandlers[0]; }

public:

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(GeantTrack *track, GeantTaskData *td);

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

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() { return "DiscreteProc"; }

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
