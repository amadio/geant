
#ifndef POSTPROPAGATIONSTAGE_H
#define POSTPROPAGATIONSTAGE_H

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
 * @brief   Special simulation stage for msc to be invoked after the transportation.
 * @class   PostPropagationStage
 * @author  M Novak
 * @date    June 2017
 *
 * Particles that have multiple scattering process will be selected at this stage.
 */

class PostPropagationStage : public Geant::SimulationStage {
public:
  /** @brief ctr */
  PostPropagationStage() {}

  /** @brief ctr */
  PostPropagationStage(Geant::GeantPropagator *prop);

  /** @brief dtr */
  ~PostPropagationStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() { return "PostPropagationStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual Geant::Handler *Select(Geant::GeantTrack *track, Geant::GeantTaskData *td);

private:

  PostPropagationStage(const PostPropagationStage &) = delete;
  PostPropagationStage &operator=(const PostPropagationStage &) = delete;

};


}     // namespace geantphysics

#endif // POSTPROPAGATIONSTAGE_H
