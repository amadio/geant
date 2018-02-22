
#ifndef PREPROPAGATIONHANDLER_H
#define PREPROPAGATIONHANDLER_H

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
 * @brief   Special handler for msc to be invoked before the transportation.
 * @class   PrePropagationHandler
 * @author  M Novak
 * @date    June 2017
 *
 * Multiple scattering process AlongStepLimitationLength method will be called, at the pre-step point (before the
 * transportation). The true physics step length limit (due to all other physics processes) will also be changed to
 * geometrical one that will be used by the transportation.
 */

class PrePropagationHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  PrePropagationHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  PrePropagationHandler(int threshold, geant::Propagator *propagator);

  /** @brief dtr */
  virtual ~PrePropagationHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::Track *track, geant::Basket& output, geant::GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket& output, geant::GeantTaskData *td);

private:
  PrePropagationHandler(const PrePropagationHandler &) = delete;
  PrePropagationHandler &operator=(const PrePropagationHandler &) = delete;

};


}      // namespace geantphysics

#endif // PREPROPAGATIONHANDLER_H
