
#ifndef POSTSTEPACTIONSTAGE_H
#define POSTSTEPACTIONSTAGE_H

// from geantV
#include "SimulationStage.h"
// from geantV
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
  class Propagator;
  class Track;
  class GeantTaskData;
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

  /** @brief Get simulation stage name */
  virtual const char *GetName() { return "PostStepActionStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::GeantTaskData *td);

private:

  PostStepActionStage(const PostStepActionStage &) = delete;
  PostStepActionStage &operator=(const PostStepActionStage &) = delete;

};

}        // namespace geantphysics

#endif   // POSTSTEPACTIONSTAGE_H
