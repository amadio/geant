
#ifndef POSTPROPAGATIONHANDLER_H
#define POSTPROPAGATIONHANDLER_H

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


class PostPropagationHandler : public Geant::Handler {
public:
  /** @brief Default constructor */
  PostPropagationHandler() : Geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  PostPropagationHandler(int threshold, Geant::GeantPropagator *propagator);

  /** @brief dtr */
  virtual ~PostPropagationHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  // virtual void DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td);

private:
  PostPropagationHandler(const PostPropagationHandler &) = delete;
  PostPropagationHandler &operator=(const PostPropagationHandler &) = delete;

};


}      // namespace geantphysics

#endif // POSTPROPAGATIONHANDLER_H
