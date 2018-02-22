
#ifndef COMPUTEINTLHANDLER_H
#define COMPUTEINTLHANDLER_H

// from geantV
#include "Handler.h"
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
 * @brief   Handler for computing the physics step limit forom normal (non msc) processes.
 * @class   ComputeIntLHandler
 * @author  M Novak
 * @date    May 2017
 */


class ComputeIntLHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  ComputeIntLHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  ComputeIntLHandler(int threshold, geant::Propagator *propagator);

  /** @brief dtr */
  virtual ~ComputeIntLHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::Track *track, geant::Basket& output, geant::TaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket& output, geant::TaskData *td);

private:
  ComputeIntLHandler(const ComputeIntLHandler &) = delete;
  ComputeIntLHandler &operator=(const ComputeIntLHandler &) = delete;

};

} // namespace geantphysics

#endif
