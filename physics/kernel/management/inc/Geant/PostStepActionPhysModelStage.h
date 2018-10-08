//===--- PostStepActionPhysModelStage.h - Geant-V ----------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PostStepActionPhysModelStage.h
 * @brief Implementation of post step action stage for vectorized EM transport
 * @author Vitalii Drohan
 */
//===----------------------------------------------------------------------===//
#ifndef POSTSTEPACTIONPHYSMODELSSTAGE_H
#define POSTSTEPACTIONPHYSMODELSSTAGE_H

// from geantV
#include "Geant/SimulationStage.h"

#include <vector>
// from geantV
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
class Propagator;
class Track;
class TaskData;
class Handler;
}
}

namespace geantphysics {

/**
 * @brief This stage filters particles based on which EM model should be sampled.
 * @class   PostStepActionPhysModelStage
 * Work is dispatched to PostStepActionPhysModelHandler
 * (based on PostStepActionStage by M. Novak)
 */

class PostStepActionPhysModelStage : public geant::SimulationStage {
public:
  /** @brief ctr */
  PostStepActionPhysModelStage() {}

  /** @brief ctr */
  PostStepActionPhysModelStage(geant::Propagator *prop);

  /** @brief dtr */
  ~PostStepActionPhysModelStage();

  VECCORE_ATT_HOST_DEVICE
  PostStepActionPhysModelStage(const PostStepActionPhysModelStage &);

  VECCORE_ATT_HOST_DEVICE
  PostStepActionPhysModelStage &operator=(const PostStepActionPhysModelStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual geant::SimulationStage *Clone() const;

  /** @brief Get simulation stage name */
  const char *GetName() const { return "PostStepActionPhysModelStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  geant::Handler *Select(geant::Track *track, geant::TaskData *td);

private:
  std::vector<geant::Handler *> fHandlersPerModel;
};

} // namespace geantphysics

#endif // POSTSTEPACTIONSTAGE_H
