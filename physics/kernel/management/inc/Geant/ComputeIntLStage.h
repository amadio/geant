
#ifndef COMPUTEINTLSTAGE_H
#define COMPUTEINTLSTAGE_H

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
  ComputeIntLStage(geant::Propagator *prop);

  /** @brief dtr */
  ~ComputeIntLStage();

  VECCORE_ATT_HOST_DEVICE
  ComputeIntLStage(const ComputeIntLStage &);

  VECCORE_ATT_HOST_DEVICE
  ComputeIntLStage &operator=(const ComputeIntLStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual geant::SimulationStage *Clone() const;

  /** @brief Get simulation stage name */
  virtual const char *GetName() const { return "ComputeIntL"; }

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  virtual int CreateHandlers();

  /** @brief Interface to select the handler matching a track */
  virtual geant::Handler *Select(geant::Track *track, geant::TaskData *td);
};

} // namespace geantphysics

#endif // COMPUTEINTLSTAGE_H
