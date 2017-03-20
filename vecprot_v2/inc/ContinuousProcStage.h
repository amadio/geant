//===--- ContinuousProcStage.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file ContinuousProcStage.h
 * @brief The stage dealing with continuous processes.
 * @details Applies continuous processes to all staged tracks
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_CONTINUOUS_PROC_STAGE
#define GEANT_CONTINUOUS_PROC_STAGE

#include "SimulationStage.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class ContinuousProcStage : public SimulationStage {

private:
  ContinuousProcStage(const ContinuousProcStage &) = delete;
  ContinuousProcStage &operator=(const ContinuousProcStage &) = delete;

public:

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(GeantTrack *track, GeantTaskData *td);

public:
  /** @brief Dummy ContinuousProcStage constructor */
  VECCORE_ATT_HOST_DEVICE
  ContinuousProcStage() {}

  /** @brief Standard ContinuousProcStage constructor */
  VECCORE_ATT_HOST_DEVICE
  ContinuousProcStage(GeantPropagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  ~ContinuousProcStage();

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() { return "ContinuousProc"; }

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
