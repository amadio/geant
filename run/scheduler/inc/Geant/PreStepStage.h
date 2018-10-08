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

#include "Geant/SimulationStage.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class PreStepStage : public SimulationStage {

public:
  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(Track *track, TaskData *td);

public:
  /** @brief Dummy PreStepStage constructor */
  VECCORE_ATT_HOST_DEVICE
  PreStepStage() {}

  /** @brief Standard PreStepStage constructor */
  VECCORE_ATT_HOST_DEVICE
  PreStepStage(Propagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~PreStepStage() {}

  VECCORE_ATT_HOST_DEVICE
  PreStepStage(const PreStepStage &);

  VECCORE_ATT_HOST_DEVICE
  PreStepStage &operator=(const PreStepStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual SimulationStage *Clone() const;

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() const { return "PreStep"; }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
