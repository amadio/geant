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
 *          stage during simulation. A stage can have one or more follow-ups. The
 *          input for a stage is a basket handled with the rask data, containing 
 *          unsorted tracks to execute the stage actions. Specialized stages have
 *          to implement at minimum the filtering selection criteria and follow-up
 *          stage after execution.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_SIMULATION_STAGE
#define GEANT_SIMULATION_STAGE

#include "Geant/Typedefs.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class GeantTrack;
class Basket;
class Filter;
class GeantPropagator;
#include "GeantFwd.h"

/** Basket processing stages. */
enum ESimulationStage {
  kUndefinedStage,       // Undefined stage type
  kSchedulingStage,      // Scheduling actions: track fetching, flushing, prioritizing
  kSampleXsecStage,      // Propose physics step by sampling total Xsec
  kGeometryStepStage,    // Compute geometry transport length
  kPropagationStage,     // Propagation in field stage
  kContinuousProcStage,  // Continuous processes stage
  kDiscreteProcStage,    // Discrete processes stage
  kBufferingStage,       // Stack-like buffering stage
  kUserStage             // Any user stage (user actions, fast simulation)
};

//template <ESimulationStage STAGE>
class SimulationStage {

  using Filters_t = vector_t<Filter *>;
  using Stages_t = vector_t<SimulationStage *>;

protected:  
  ESimulationStage fType = kUndefinedStage; ///< Processing stage type
  GeantPropagator *fPropagator = nullptr;   ///< Propagator owning this stage
  int fId = -1;                             ///< Unique stage id
  Filters_t fFilters;                       ///< Array of filters
  Stages_t fFollowUps;                      ///< Follow-up stages for processed tracks
  
private:
  SimulationStage(const SimulationStage &) = delete;
  SimulationStage &operator=(const SimulationStage &) = delete;

// The functions below are the interfaces for derived simulation stages.
protected:

  /** @brief Interface to create all filters for the simulation stage
   *  @return Number of filters created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateFilters() { return 0; }

public:

  /** @brief Interface to select the filter matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Filter *SelectFilter(GeantTrack *track) = 0;

  /** @brief Interface to select the next stage for a processed track, from the list of follow-ups */
  VECCORE_ATT_HOST_DEVICE
  virtual int SelectFollowUp(GeantTrack *) { return 0; }
  
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

//=== The stage processing methods === //

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

  /** @brief Getter for type */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  ESimulationStage GetType() const { return fType; }

  /** @brief Getter for type */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int NFollowUps() const { return fFollowUps.size(); }

  /** @brief Add next filter */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void AddFilter(Filter *filter);

  /** @brief Getter for number of filters */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNfilters() const { return fFilters.size(); }

  /** @brief Add a follow-up stage */
  VECCORE_ATT_HOST_DEVICE
  void AddFollowUpStage(SimulationStage *stage)
  {
    // Not allowed that the same stage is a follow-up
    assert(stage != this);
    fFollowUps.push_back(stage);
  }

  /** @brief Getter for type */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  SimulationStage *GetFollowUpStage(int istage = 0) const { return fFollowUps[istage]; }

  /** @brief Fork output basket to smaller baskets if stage has more daughters */
  VECCORE_ATT_HOST_DEVICE
  void ForkToFollowUps(Basket *output, GeantTaskData *td);
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
