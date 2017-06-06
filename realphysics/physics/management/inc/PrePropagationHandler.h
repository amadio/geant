
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
  // virtual void DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td);

private:
  PrePropagationHandler(const PrePropagationHandler &) = delete;
  PrePropagationHandler &operator=(const PrePropagationHandler &) = delete;

};


}      // namespace geantphysics

#endif // PREPROPAGATIONHANDLER_H
