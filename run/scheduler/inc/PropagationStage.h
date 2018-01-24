//===--- PropagationStage.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PropagationStage.h
 * @brief The track propagaion as simulation stage.
 * @details This simulation stage deals with propagating tracks in a single volume
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PROPAGATION_STAGE
#define GEANT_PROPAGATION_STAGE

#include "SimulationStage.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class PropagationStage : public SimulationStage {

protected:
  bool fHasField = true;         ///< Setup has field
private:
  PropagationStage(const PropagationStage &) = delete;
  PropagationStage &operator=(const PropagationStage &) = delete;

public:
  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(GeantTrack *track, GeantTaskData *td);

public:
  /** @brief Dummy PropagationStage constructor */
  VECCORE_ATT_HOST_DEVICE
  PropagationStage() {}

  /** @brief Standard PropagationStage constructor */
  VECCORE_ATT_HOST_DEVICE
  PropagationStage(GeantPropagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  ~PropagationStage();

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() { return "Propagation"; }

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
