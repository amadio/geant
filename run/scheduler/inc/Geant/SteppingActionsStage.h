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

#include "Geant/SimulationStage.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class SteppingActionsStage : public SimulationStage {

public:
  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(Track *track, TaskData *td);

public:
  /** @brief Dummy SteppingActionsStage constructor */
  VECCORE_ATT_HOST_DEVICE
  SteppingActionsStage() {}

  /** @brief Standard SteppingActionsStage constructor */
  VECCORE_ATT_HOST_DEVICE
  SteppingActionsStage(Propagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~SteppingActionsStage() {}

  VECCORE_ATT_HOST_DEVICE
  SteppingActionsStage(const SteppingActionsStage &);

  VECCORE_ATT_HOST_DEVICE
  SteppingActionsStage &operator=(const SteppingActionsStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual SimulationStage *Clone() const;

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() const { return "SteppingActions"; }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
