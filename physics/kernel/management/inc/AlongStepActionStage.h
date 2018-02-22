
#ifndef ALONGSTEPACTIONSTAGE_H
#define ALONGSTEPACTIONSTAGE_H

// from geantV
#include "SimulationStage.h"
// from geantV
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
  class GeantPropagator;
  class GeantTrack;
  class GeantTaskData;
  class Handler;
}
}


namespace geantphysics {

/**
 * @brief   Simulation stage for along-step-actions.
 * @class   AlongStepActionStage
 * @author  M Novak
 * @date    May 2017
 */


class AlongStepActionStage : public geant::SimulationStage {
public:

  /** @brief ctr */
  AlongStepActionStage() {}

  /** @brief ctr */
  AlongStepActionStage(geant::GeantPropagator *prop);

  /** @brief dtr */
  ~AlongStepActionStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() { return "AlongStepActionStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::GeantTrack *track, geant::GeantTaskData *td);

private:

  AlongStepActionStage(const AlongStepActionStage &) = delete;
  AlongStepActionStage &operator=(const AlongStepActionStage &) = delete;

};

}       // namespace geantphysics

#endif  // ALONGSTEPACTIONSTAGE_H
