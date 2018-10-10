
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
} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant

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

  VECCORE_ATT_HOST_DEVICE
  FastSimStage(const FastSimStage &);

  VECCORE_ATT_HOST_DEVICE
  FastSimStage &operator=(const FastSimStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual geant::SimulationStage *Clone() const;

  /** @brief Get simulation stage name */
  virtual const char *GetName() const { return "FastSim"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::TaskData *td);
};

} // namespace geantphysics

#endif // FASTSIMSTAGE_H
