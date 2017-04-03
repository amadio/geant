//===--- SimulationStage.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file SimulationStage.h
 * @brief Abstraction for a simulation stage.
 * @details A simulation stage contains all the utilities allowing selecting,
 *          basketizing and running specific actions for a given track processing
 *          stage during simulation. A stage can have one or more follow-ups. The
 *          input for a stage is a basket handled with the rask data, containing 
 *          unsorted tracks to execute the stage actions. Specialized stages have
 *          to implement at minimum the selection criteria and follow-up
 *          stage after execution.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_SIMULATION_STAGE
#define GEANT_SIMULATION_STAGE

#include "Geant/Typedefs.h"
#include "GeantTrack.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class Basket;
class Handler;
class GeantPropagator;
#include "GeantFwd.h"

//template <ESimulationStage STAGE>
class SimulationStage {

  using Handlers_t = vector_t<Handler *>;

protected:  
  ESimulationStage fType = kPreStepStage;   ///< Processing stage type
  GeantPropagator *fPropagator = nullptr;   ///< Propagator owning this stage
  int fId = 0;                              ///< Unique stage id
  int fUserActionsStage = 0;                ///< User actions stage to be executed right after
  int fFollowUpStage = 0;                   ///< In case there is a single follow-up store its id
  bool fEndStage = false;                   ///< Marker for stage at end of stepping

  Handlers_t fHandlers;                   ///< Array of handlers for the stage
  
private:
  SimulationStage(const SimulationStage &) = delete;
  SimulationStage &operator=(const SimulationStage &) = delete;

  VECCORE_ATT_HOST_DEVICE
  int CopyToFollowUps(Basket &output, GeantTaskData *td);
// The functions below are the interfaces for derived simulation stages.

public:

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers() = 0;

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(GeantTrack *track, GeantTaskData *td) = 0;

public:
  /** @brief Dummy SimulationStage constructor */
  VECCORE_ATT_HOST_DEVICE
  SimulationStage() {}

  /** @brief Standard SimulationStage constructor */
  VECCORE_ATT_HOST_DEVICE
  SimulationStage(ESimulationStage type, GeantPropagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  ~SimulationStage();

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() { return nullptr; }
  

//=== The stage processing methods === //

  /** @brief Process a basket of tracks marked for the stage
   *  @return Number of tracks processed
   */
  VECCORE_ATT_HOST_DEVICE
  virtual void ActivateBasketizing(bool) {}

  /** @brief Process a basket of tracks marked for the stage
   *  @return Number of tracks processed
   */
  VECCORE_ATT_HOST_DEVICE
  int Process(GeantTaskData *td);

  /** @brief Flush all tracks from the simulation stage basketizers and execute stage
   *  @return Number of tracks flushed
   */
  VECCORE_ATT_HOST_DEVICE
  int FlushAndProcess(GeantTaskData *td);

  /** @brief Getter for type */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  ESimulationStage GetType() const { return fType; }

  /** @brief Getter for Id */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetId() const { return fId; }

  /** @brief Set follow-up stage */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetFollowUpStage(ESimulationStage stage) { fFollowUpStage = (int)stage; } 

  /** @brief Getter for type */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetFollowUpStage() const { return fFollowUpStage; }

  /** @brief Add next handler */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void AddHandler(Handler *handler) { fHandlers.push_back(handler); }

  /** @brief Getter for number of handlers */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNhandlers() const { return fHandlers.size(); }

  /** @brief Getter for a given handler */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Handler *GetHandler(int i) const { return fHandlers[i]; }

  /** @brief Set user actions stage */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetUserActionsStage(ESimulationStage stage) { fUserActionsStage = (int)stage; } 

  /** @brief Getter for user actions stage */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetUserActionsStage() const { return fUserActionsStage; } 

  /** @brief Setter for end stage */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetEndStage() { fEndStage = true; } 

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
