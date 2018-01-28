
#ifndef ATRESTACTIONHANDLER_H
#define ATRESTACTIONHANDLER_H

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
 * @brief   Handler for at-rest actions.
 * @class   AtRestActionHandler
 * @author  M Novak
 * @date    January 2018
 */

class AtRestActionHandler : public Geant::Handler {
public:
  /** @brief Default constructor */
  AtRestActionHandler() : Geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  AtRestActionHandler(int threshold, Geant::GeantPropagator *propagator);

  /** @brief dtr */
  virtual ~AtRestActionHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(Geant::GeantTrack *track, Geant::Basket& output, Geant::GeantTaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(Geant::Basket &input, Geant::Basket& output, Geant::GeantTaskData *td);

private:
  AtRestActionHandler(const AtRestActionHandler &) = delete;
  AtRestActionHandler &operator=(const AtRestActionHandler &) = delete;

};

}        // namespace geantphysics

#endif   // ATRESTACTIONHANDLER_H
