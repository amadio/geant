//===--- PreStepStage.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PreStepStage.h
 * @brief The pre-step stage triggers user pre-step actions.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PRE_STEP_STAGE
#define GEANT_PRE_STEP_STAGE

#include "SimulationStage.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class PreStepStage : public SimulationStage {

private:
  PreStepStage(const PreStepStage &) = delete;
  PreStepStage &operator=(const PreStepStage &) = delete;

public:
  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(GeantTrack *track, GeantTaskData *td);

public:
  /** @brief Dummy PreStepStage constructor */
  VECCORE_ATT_HOST_DEVICE
  PreStepStage() {}

  /** @brief Standard PreStepStage constructor */
  VECCORE_ATT_HOST_DEVICE
  PreStepStage(GeantPropagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  ~PreStepStage();

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() { return "PreStep"; }

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
