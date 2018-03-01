
#ifndef POSTSTEPACTIONPHYSPROCESSHANDLER_H
#define POSTSTEPACTIONPHYSPROCESSHANDLER_H

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

/**
 * @brief   Handler for post-step actions specific for physics process (discrete interactions).
 * @class   PostStepActionPhysProcessHandler
 */

class PostStepActionPhysProcessHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  PostStepActionPhysProcessHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  PostStepActionPhysProcessHandler(int threshold, geant::Propagator *propagator);

  /** @brief dtr */
  virtual ~PostStepActionPhysProcessHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td);

private:
  PostStepActionPhysProcessHandler(const PostStepActionPhysProcessHandler &) = delete;
  PostStepActionPhysProcessHandler &operator=(const PostStepActionPhysProcessHandler &) = delete;
};

} // namespace geantphysics

#endif // POSTSTEPACTIONHANDLER_H
