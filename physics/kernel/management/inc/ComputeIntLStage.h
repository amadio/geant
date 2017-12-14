
#ifndef COMPUTEINTLSTAGE_H
#define COMPUTEINTLSTAGE_H

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
 * @brief   Simulation stage for normal (non msc) physics step limit computation.
 * @class   ComputeIntLStage
 * @author  M Novak
 * @date    May 2017
 */


class ComputeIntLStage : public Geant::SimulationStage {
public:
  /** @brief ctr */
  ComputeIntLStage() {}

  /** @brief ctr */
  ComputeIntLStage(Geant::GeantPropagator *prop);

  /** @brief dtr */
  ~ComputeIntLStage();

  /** @brief Get simulation stage name */
  virtual const char *GetName() { return "ComputeIntLStage"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual Geant::Handler *Select(Geant::GeantTrack *track, Geant::GeantTaskData *td);


private:

  ComputeIntLStage(const ComputeIntLStage &) = delete;
  ComputeIntLStage &operator=(const ComputeIntLStage &) = delete;

// ?
//  Handler *GetHandler(int) { return fHandlers[0]; }
};

}      // namespace geantphysics

#endif // COMPUTEINTLSTAGE_H
