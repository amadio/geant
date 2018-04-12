
#ifndef POSTPROPAGATIONSTAGE_H
#define POSTPROPAGATIONSTAGE_H

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
 * @brief   Special simulation stage for msc to be invoked after the transportation.
 * @class   PostPropagationStage
 * @author  M Novak
 * @date    June 2017
 *
 * Particles that have multiple scattering process will be selected at this stage.
 */

class PostPropagationStage : public geant::SimulationStage {
public:
  /** @brief ctr */
  PostPropagationStage() {}

  /** @brief ctr */
  PostPropagationStage(geant::Propagator *prop);

  /** @brief dtr */
  ~PostPropagationStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() const { return "PostPropagation"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::TaskData *td);

private:
  PostPropagationStage(const PostPropagationStage &) = delete;
  PostPropagationStage &operator=(const PostPropagationStage &) = delete;
};

} // namespace geantphysics

#endif // POSTPROPAGATIONSTAGE_H
