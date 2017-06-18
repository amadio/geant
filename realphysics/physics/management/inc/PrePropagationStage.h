
#ifndef PREPROPAGATIONSTAGE_H
#define PREPROPAGATIONSTAGE_H

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
 * @brief   Special simulation stage for msc to be invoked before the transportation.
 * @class   PrePropagationStage
 * @author  M Novak
 * @date    June 2017
 *
 * Particles that have multiple scattering process will be selected at this stage.
 */


class PrePropagationStage : public Geant::SimulationStage {
public:
  /** @brief ctr */
  PrePropagationStage() {}

  /** @brief ctr */
  PrePropagationStage(Geant::GeantPropagator *prop);

  /** @brief dtr */
  ~PrePropagationStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() { return "PrePropagationStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual Geant::Handler *Select(Geant::GeantTrack *track, Geant::GeantTaskData *td);

private:

  PrePropagationStage(const PrePropagationStage &) = delete;
  PrePropagationStage &operator=(const PrePropagationStage &) = delete;

};


}     // namespace geantphysics

#endif // PREPROPAGATIONSTAGE_H
