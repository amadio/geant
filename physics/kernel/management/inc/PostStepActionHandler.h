
#ifndef POSTSTEPACTIONHANDLER_H
#define POSTSTEPACTIONHANDLER_H

// from geantV
#include "Handler.h"
// from geantV
namespace geant {
  inline namespace GEANT_IMPL_NAMESPACE {
  class GeantPropagator;
  class GeantTrack;
  class GeantTaskData;
  class Basket;
}
}


namespace geantphysics {

/**
 * @brief   Handler for post-step actions (discrete interactions).
 * @class   PostStepActionHandler
 * @author  M Novak
 * @date    May 2017
 */

class PostStepActionHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  PostStepActionHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  PostStepActionHandler(int threshold, geant::GeantPropagator *propagator);

  /** @brief dtr */
  virtual ~PostStepActionHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::GeantTrack *track, geant::Basket& output, geant::GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket& output, geant::GeantTaskData *td);

private:
  PostStepActionHandler(const PostStepActionHandler &) = delete;
  PostStepActionHandler &operator=(const PostStepActionHandler &) = delete;

};

}        // namespace geantphysics

#endif   // POSTSTEPACTIONHANDLER_H
