
#ifndef POSTSTEPACTIONPHYSPROCESSSTAGE_H
#define POSTSTEPACTIONPHYSPROCESSSTAGE_H

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
 * @brief   Simulation stage to select particles that post-step action based on their physics process (discrete
 * interaction) need to be invoked.
 * @class   PostStepActionPhysProcessStage
 */

class PostStepActionPhysModelStage : public geant::SimulationStage {
public:
  /** @brief ctr */
  PostStepActionPhysModelStage() {}

  /** @brief ctr */
  PostStepActionPhysModelStage(geant::Propagator *prop);

  /** @brief dtr */
  ~PostStepActionPhysModelStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() { return "PostStepActionPhysModelStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::TaskData *td);

private:
  PostStepActionPhysModelStage(const PostStepActionPhysModelStage &) = delete;
  PostStepActionPhysModelStage &operator=(const PostStepActionPhysModelStage &) = delete;

  std::vector<geant::Handler *> fHandlersPerModel;
};

} // namespace geantphysics

#endif // POSTSTEPACTIONSTAGE_H
