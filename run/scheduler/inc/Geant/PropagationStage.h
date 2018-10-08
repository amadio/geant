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

#include "Geant/SimulationStage.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class PropagationStage : public SimulationStage {

protected:
  bool fHasField = true;           ///< Setup has field
  bool fLocalFieldHandler = false; ///< Field handler is thread local

public:
  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(Track *track, TaskData *td);

public:
  /** @brief Dummy PropagationStage constructor */
  VECCORE_ATT_HOST_DEVICE
  PropagationStage() {}

  /** @brief Standard PropagationStage constructor */
  VECCORE_ATT_HOST_DEVICE
  PropagationStage(Propagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~PropagationStage() {}

  VECCORE_ATT_HOST_DEVICE
  PropagationStage(const PropagationStage &);

  VECCORE_ATT_HOST_DEVICE
  PropagationStage &operator=(const PropagationStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual SimulationStage *Clone() const;

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() const { return "Propagation"; }

  /** @brief Setter for local field handler */
  VECCORE_ATT_HOST_DEVICE
  void SetLocalFieldHandler(bool flag = true) { fLocalFieldHandler = flag; }

  /** @brief Getter for local field handler */
  VECCORE_ATT_HOST_DEVICE
  bool GetLocalFieldHandler() const { return fLocalFieldHandler; }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
