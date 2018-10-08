
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

  VECCORE_ATT_HOST_DEVICE
  PostPropagationStage(const PostPropagationStage &);

  VECCORE_ATT_HOST_DEVICE
  PostPropagationStage &operator=(const PostPropagationStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual geant::SimulationStage *Clone() const;

  /** @brief Get simulation stage name */
  virtual const char *GetName() const { return "PostPropagation"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::TaskData *td);

};

} // namespace geantphysics

#endif // POSTPROPAGATIONSTAGE_H
