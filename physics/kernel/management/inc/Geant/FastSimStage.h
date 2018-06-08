
#ifndef FASTSIMSTAGE_H
#define FASTSIMSTAGE_H

// from geantV
#include "Geant/SimulationStage.h"
// from geantV
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
class Propagator;
class Track;
class TaskData;
class Handler;
}
}

namespace geantphysics {

/**
 * @brief   Simulation stage Fast Simulation process.
 * @class   FastSimStage
 * @author  W Pokorski
 * @date    May 2018
 */

class FastSimStage : public geant::SimulationStage {
public:
  /** @brief ctr */
  FastSimStage() {}

  /** @brief ctr */
  FastSimStage(geant::Propagator *prop);

  /** @brief dtr */
  ~FastSimStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() const { return "FastSim"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::TaskData *td);

private:
  FastSimStage(const FastSimStage &) = delete;
  FastSimStage &operator=(const FastSimStage &) = delete;

  // ?
  //  Handler *GetHandler(int) { return fHandlers[0]; }
};

}// namespace geantphysics

#endif // FASTSIMSTAGE_H
