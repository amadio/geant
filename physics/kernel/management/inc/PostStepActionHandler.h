
#ifndef POSTSTEPACTIONHANDLER_H
#define POSTSTEPACTIONHANDLER_H

// from geantV
#include "Handler.h"
// from geantV
namespace Geant {
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

class PostStepActionHandler : public Geant::Handler {
public:
  /** @brief Default constructor */
  PostStepActionHandler() : Geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  PostStepActionHandler(int threshold, Geant::GeantPropagator *propagator);

  /** @brief dtr */
  virtual ~PostStepActionHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td);

private:
  PostStepActionHandler(const PostStepActionHandler &) = delete;
  PostStepActionHandler &operator=(const PostStepActionHandler &) = delete;

};

}        // namespace geantphysics

#endif   // POSTSTEPACTIONHANDLER_H
