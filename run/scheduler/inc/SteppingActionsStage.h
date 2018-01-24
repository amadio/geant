//===--- SteppingActionsStage.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file SteppingActionsStage.h
 * @brief Stepping actions as simulation stage.
 * @details This simulation stage deals with invoking the user stepping actions.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_STEPPING_ACTIONS_STAGE
#define GEANT_STEPPING_ACTIONS_STAGE

#include "SimulationStage.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class SteppingActionsStage : public SimulationStage {

private:
  SteppingActionsStage(const SteppingActionsStage &) = delete;
  SteppingActionsStage &operator=(const SteppingActionsStage &) = delete;

public:
  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(GeantTrack *track, GeantTaskData *td);

public:
  /** @brief Dummy SteppingActionsStage constructor */
  VECCORE_ATT_HOST_DEVICE
  SteppingActionsStage() {}

  /** @brief Standard SteppingActionsStage constructor */
  VECCORE_ATT_HOST_DEVICE
  SteppingActionsStage(GeantPropagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  ~SteppingActionsStage();

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() { return "SteppingActions"; }

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
