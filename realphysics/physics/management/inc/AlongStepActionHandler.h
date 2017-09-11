
#ifndef ALONGSTEPACTIONHANDLER_H
#define ALONGSTEPACTIONHANDLER_H

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
 * @brief   Handler for along-step-actions.
 * @class   AlongStepActionHandler
 * @author  M Novak
 * @date    May 2017
 */

class AlongStepActionHandler : public Geant::Handler {
public:
  /** @brief Default constructor */
  AlongStepActionHandler() : Geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  AlongStepActionHandler(int threshold, Geant::GeantPropagator *propagator);

  /** @brief dtr */
  virtual ~AlongStepActionHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td);

private:
  AlongStepActionHandler(const AlongStepActionHandler &) = delete;
  AlongStepActionHandler &operator=(const AlongStepActionHandler &) = delete;

};

}        // namespace geantphysics

#endif   // ALONGSTEPACTIONHANDLER_H
