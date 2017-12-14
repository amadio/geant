
#ifndef COMPUTEINTLHANDLER_H
#define COMPUTEINTLHANDLER_H

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
 * @brief   Handler for computing the physics step limit forom normal (non msc) processes.
 * @class   ComputeIntLHandler
 * @author  M Novak
 * @date    May 2017
 */


class ComputeIntLHandler : public Geant::Handler {
public:
  /** @brief Default constructor */
  ComputeIntLHandler() : Geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  ComputeIntLHandler(int threshold, Geant::GeantPropagator *propagator);

  /** @brief dtr */
  virtual ~ComputeIntLHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td);

private:
  ComputeIntLHandler(const ComputeIntLHandler &) = delete;
  ComputeIntLHandler &operator=(const ComputeIntLHandler &) = delete;

};

} // namespace geantphysics

#endif
