
#ifndef FASTSIMHANDLER_H
#define FASTSIMHANDLER_H

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
 * @brief   Handler for fast simulation processes.
 * @class   FastSimHandler
 * @author  W Pokorski
 * @date    May 2018
 */

class FastSimHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  FastSimHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  FastSimHandler(int threshold, geant::Propagator *propagator);

  /** @brief dtr */
  virtual ~FastSimHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td);

private:
  FastSimHandler(const FastSimHandler &) = delete;
  FastSimHandler &operator=(const FastSimHandler &) = delete;
};

} // namespace geantphysics

#endif
