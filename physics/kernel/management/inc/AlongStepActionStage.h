
#ifndef ALONGSTEPACTIONSTAGE_H
#define ALONGSTEPACTIONSTAGE_H

// from geantV
#include "SimulationStage.h"
// from geantV
namespace Geant {
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


class AlongStepActionStage : public Geant::SimulationStage {
public:

  /** @brief ctr */
  AlongStepActionStage() {}

  /** @brief ctr */
  AlongStepActionStage(Geant::GeantPropagator *prop);

  /** @brief dtr */
  ~AlongStepActionStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() { return "AlongStepActionStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual Geant::Handler *Select(Geant::GeantTrack *track, Geant::GeantTaskData *td);

private:

  AlongStepActionStage(const AlongStepActionStage &) = delete;
  AlongStepActionStage &operator=(const AlongStepActionStage &) = delete;

};

}       // namespace geantphysics

#endif  // ALONGSTEPACTIONSTAGE_H
