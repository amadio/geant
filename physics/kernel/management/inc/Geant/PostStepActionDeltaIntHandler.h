#ifndef POSTSTEPACTIONDELTAINTHANDLER_H
#define POSTSTEPACTIONDELTAINTHANDLER_H

// from geantV
#include "Geant/Handler.h"
// from geantV
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
class Propagator;
class Track;
class TaskData;
class Basket;
}
}

namespace geantphysics {
class EMModel;
/**
 * @brief   Handler for post-step actions specific for physics process (discrete interactions).
 * @class   PostStepActionDeltaIntHandler
 */

class PostStepActionDeltaIntHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  PostStepActionDeltaIntHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  PostStepActionDeltaIntHandler(int threshold, geant::Propagator *propagator);

  /** @brief dtr */
  virtual ~PostStepActionDeltaIntHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td);

private:
  PostStepActionDeltaIntHandler(const PostStepActionDeltaIntHandler &) = delete;
  PostStepActionDeltaIntHandler &operator=(const PostStepActionDeltaIntHandler &) = delete;
};

} // namespace geantphysics

#endif // POSTSTEPACTIONHANDLER_H
