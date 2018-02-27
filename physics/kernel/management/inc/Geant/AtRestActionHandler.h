
#ifndef ATRESTACTIONHANDLER_H
#define ATRESTACTIONHANDLER_H

// from geantV
#include "Geant/Handler.h"
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
 * @brief   Handler for at-rest actions.
 * @class   AtRestActionHandler
 * @author  M Novak
 * @date    January 2018
 */

class AtRestActionHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  AtRestActionHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  AtRestActionHandler(int threshold, geant::Propagator *propagator);

  /** @brief dtr */
  virtual ~AtRestActionHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::Track *track, geant::Basket& output, geant::TaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket& output, geant::TaskData *td);

private:
  AtRestActionHandler(const AtRestActionHandler &) = delete;
  AtRestActionHandler &operator=(const AtRestActionHandler &) = delete;

};

}        // namespace geantphysics

#endif   // ATRESTACTIONHANDLER_H
