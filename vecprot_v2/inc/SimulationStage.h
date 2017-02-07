//===--- SimulationStage.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file SimulationStage.h
 * @brief Abstraction for a simulation stage.
 * @details A simulation stage contains all the utilities allowing filtering,
 *          basketizing and running specific actions for a given track processing
 *          stage during simulation. Templated class allowing specializing for
 *          different stage types, such as scheduling, x-sec sampling, geometry,
 *          propagation or physics along/post-stepping actions.
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_SIMULATION_STAGE
#define GEANT_SIMULATION_STAGE

#include "Geant/Typedefs.h"
#include "priority_queue.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class GeantTrack;
class Basket;
class Filter;
#include "GeantFwd.h"

/** Basket processing stages. */
enum EStage {
  kScheduling,      // Scheduling actions: track fetching, flushing, prioritizing
  kSampleXsec,      // Propose physics step by sampling total Xsec
  kTransportLength, // Compute geometry transport length
  kPropagate,       // Propagation in field
  kAlongStep,       // Along step actions (e.g. continuous energy loss)
  kPostStep,        // Post step actions
  kStackLikeBuffer, // Stack-like buffering
  kUserActions      // User actions
};

//template <EStage STAGE>
class SimulationStage {

  using queue_t = priority_queue<Basket *>;
protected:  
  EStage fType = kScheduling;          ///< Locality type
  int fNfilters = 0;                   ///< Number of filters
  int fCapacity = 128;                 ///< Size of the array of filters
  Filter **fFilters = nullptr;         ///< Array of filters
  queue_t *fToProcess = nullptr;       ///< Queue of baskets to process
  SimulationStage *fNextStage = nullptr; ///< Next simulation stage
  
private:
  SimulationStage(const SimulationStage &) = delete;
  SimulationStage &operator=(const SimulationStage &) = delete;

protected:
  /** @brief Interface to create all filters for the simulation stage
   *  @return Number of filters created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateFilters() { return 0; }

public:
  /** @brief Interface to select the filter matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Filter *Select(GeantTrack *track) = 0;

  
public:
  /** @brief Default SimulationStage constructor */
  VECCORE_ATT_HOST_DEVICE
  SimulationStage() {}

  /** @brief SimulationStage constructor for a given type */
  VECCORE_ATT_HOST_DEVICE
  SimulationStage(EStage type);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  ~SimulationStage();

  /** @brief Getter for type */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  EStage GetType() const { return fType; }

  /** @brief Add next filter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void AddFilter(Filter *filter);

  /** @brief Add a basket to the input queue */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void AddBasket(Basket *btodo) { fToProcess->push(btodo); }

  /** @brief Getter for the feeder queue */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  queue_t *GetInputQ() const { return fToProcess; }

  /** @brief Process a basket of tracks marked for the stage
   *  @return Number of tracks processed
   */
  VECCORE_ATT_HOST_DEVICE
  int Process(Basket &input, GeantTaskData *td);

  /** @brief Flush all tracks from the simulation stage basketizers and execute stage
   *  @return Number of tracks flushed
   */
  VECCORE_ATT_HOST_DEVICE
  int FlushAndProcess(Basket &btodo, GeantTaskData *td);

  /** @brief Setter for next stage to be completed after this one */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetNextStage(SimulationStage *stage) { fNextStage = stage; }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
