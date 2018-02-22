
#ifndef COMPUTEINTLSTAGE_H
#define COMPUTEINTLSTAGE_H

// from geantV
#include "SimulationStage.h"
// from geantV
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
  class GeantPropagator;
  class GeantTrack;
  class GeantTaskData;
  class Handler;
}
}


namespace geantphysics {

/**
 * @brief   Simulation stage for normal (non msc) physics step limit computation.
 * @class   ComputeIntLStage
 * @author  M Novak
 * @date    May 2017
 */


class ComputeIntLStage : public geant::SimulationStage {
public:
  /** @brief ctr */
  ComputeIntLStage() {}

  /** @brief ctr */
  ComputeIntLStage(geant::GeantPropagator *prop);

  /** @brief dtr */
  ~ComputeIntLStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() { return "ComputeIntLStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::GeantTrack *track, geant::GeantTaskData *td);


private:

  ComputeIntLStage(const ComputeIntLStage &) = delete;
  ComputeIntLStage &operator=(const ComputeIntLStage &) = delete;

// ?
//  Handler *GetHandler(int) { return fHandlers[0]; }
};

}      // namespace geantphysics

#endif // COMPUTEINTLSTAGE_H
