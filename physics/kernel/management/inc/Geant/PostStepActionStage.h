
#ifndef POSTSTEPACTIONSTAGE_H
#define POSTSTEPACTIONSTAGE_H

// from geantV
#include "Geant/SimulationStage.h"
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
 * @brief   Simulation stage to select particles that post-step action (discrete interaction) need to be invoked.
 * @class   PostStepActionStage
 * @author  M Novak
 * @date    May 2017
 */

class PostStepActionStage : public geant::SimulationStage {
public:
  /** @brief ctr */
  PostStepActionStage() {}

  /** @brief ctr */
  PostStepActionStage(geant::Propagator *prop);

  /** @brief dtr */
  ~PostStepActionStage();

  VECCORE_ATT_HOST_DEVICE
  PostStepActionStage(const PostStepActionStage &);

  VECCORE_ATT_HOST_DEVICE
  PostStepActionStage &operator=(const PostStepActionStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual geant::SimulationStage *Clone() const;

  /** @brief Get simulation stage name */
  virtual const char *GetName() const { return "PostStepAction"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::TaskData *td);

};

} // namespace geantphysics

#endif // POSTSTEPACTIONSTAGE_H
