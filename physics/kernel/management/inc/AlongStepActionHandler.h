
#ifndef ALONGSTEPACTIONHANDLER_H
#define ALONGSTEPACTIONHANDLER_H

// from geantV
#include "Handler.h"
// from geantV
namespace geant {
  inline namespace GEANT_IMPL_NAMESPACE {
  class Propagator;
  class Track;
  class GeantTaskData;
  class Basket;
}
}


namespace geantphysics {

/**
 * @brief   Handler for along-step-actions.
 * @class   AlongStepActionHandler
 * @author  M Novak
 * @date    May 2017
 */

class AlongStepActionHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  AlongStepActionHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  AlongStepActionHandler(int threshold, geant::Propagator *propagator);

  /** @brief dtr */
  virtual ~AlongStepActionHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::Track *track, geant::Basket& output, geant::GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket& output, geant::GeantTaskData *td);

private:
  AlongStepActionHandler(const AlongStepActionHandler &) = delete;
  AlongStepActionHandler &operator=(const AlongStepActionHandler &) = delete;

};

}        // namespace geantphysics

#endif   // ALONGSTEPACTIONHANDLER_H
