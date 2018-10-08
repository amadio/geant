
#ifndef ATRESTACTIONSTAGE_H
#define ATRESTACTIONSTAGE_H

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
 * @brief   Simulation stage to select particles that at-rest interaction needs to be invoked.
 * @class   AtRestActionStage
 * @author  M Novak
 * @date    January 2018
 */

class AtRestActionStage : public geant::SimulationStage {
public:
  /** @brief ctr */
  AtRestActionStage() {}

  /** @brief ctr */
  AtRestActionStage(geant::Propagator *prop);

  /** @brief dtr */
  ~AtRestActionStage();

  VECCORE_ATT_HOST_DEVICE
  AtRestActionStage(const AtRestActionStage &);

  VECCORE_ATT_HOST_DEVICE
  AtRestActionStage &operator=(const AtRestActionStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual geant::SimulationStage *Clone() const;

  /** @brief Get simulation stage name */
  virtual const char *GetName() const { return "AtRestAction"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::TaskData *td);

};

} // namespace geantphysics

#endif // ATRESTACTIONSTAGE_H
