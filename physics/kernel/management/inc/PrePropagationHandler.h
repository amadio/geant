
#ifndef PREPROPAGATIONHANDLER_H
#define PREPROPAGATIONHANDLER_H

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
 * @brief   Special handler for msc to be invoked before the transportation.
 * @class   PrePropagationHandler
 * @author  M Novak
 * @date    June 2017
 *
 * Multiple scattering process AlongStepLimitationLength method will be called, at the pre-step point (before the
 * transportation). The true physics step length limit (due to all other physics processes) will also be changed to
 * geometrical one that will be used by the transportation.
 */

class PrePropagationHandler : public Geant::Handler {
public:
  /** @brief Default constructor */
  PrePropagationHandler() : Geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  PrePropagationHandler(int threshold, Geant::GeantPropagator *propagator);

  /** @brief dtr */
  virtual ~PrePropagationHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td);

private:
  PrePropagationHandler(const PrePropagationHandler &) = delete;
  PrePropagationHandler &operator=(const PrePropagationHandler &) = delete;

};


}      // namespace geantphysics

#endif // PREPROPAGATIONHANDLER_H
