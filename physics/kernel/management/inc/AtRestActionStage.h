
#ifndef ATRESTACTIONSTAGE_H
#define ATRESTACTIONSTAGE_H

// from geantV
#include "SimulationStage.h"
// from geantV
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
  class GeantPropagator;
  class GeantTrack;
  class GeantTaskData;
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

class AtRestActionStage : public Geant::SimulationStage {
public:

  /** @brief ctr */
  AtRestActionStage() {}

  /** @brief ctr */
  AtRestActionStage(Geant::GeantPropagator *prop);

  /** @brief dtr */
  ~AtRestActionStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() { return "AtRestActionStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual Geant::Handler *Select(Geant::GeantTrack *track, Geant::GeantTaskData *td);

private:

  AtRestActionStage(const AtRestActionStage &) = delete;
  AtRestActionStage &operator=(const AtRestActionStage &) = delete;

};

}        // namespace geantphysics

#endif   // ATRESTACTIONSTAGE_H
