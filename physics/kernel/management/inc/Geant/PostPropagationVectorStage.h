
#ifndef POSTPROPAGATIONVECTORSTAGE_H
#define POSTPROPAGATIONVECTORSTAGE_H

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

class PostPropagationVectorStage : public geant::SimulationStage {
public:
  /** @brief ctr */
  PostPropagationVectorStage() {}

  /** @brief ctr */
  PostPropagationVectorStage(geant::Propagator *prop);

  /** @brief dtr */
  ~PostPropagationVectorStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() const { return "PostPropagationVector"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::TaskData *td);

private:
  PostPropagationVectorStage(const PostPropagationVectorStage &) = delete;
  PostPropagationVectorStage &operator=(const PostPropagationVectorStage &) = delete;

  std::vector<geant::Handler *> fHandlersPerModel;
};

} // namespace geantphysics

#endif // POSTPROPAGATIONSTAGE_H
