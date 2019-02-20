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
#include "Geant/Track.h"
#include "Geant/Handler.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class TaskData;
class Basket;
class Propagator;
#include "Geant/Fwd.h"

// template <ESimulationStage STAGE>
class SimulationStage {

#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  // On cuda there is one propagator per thread.  So (for now), no need
  // for atomics.
  template <typename T>
  using atomic_t = T;
#else
  template <typename T>
  using atomic_t = std::atomic<T>;
#endif
  using Handlers_t = vector_t<Handler *>;

protected:
  ESimulationStage fType  = kPreStepStage; ///< Processing stage type
  Propagator *fPropagator = nullptr;       ///< Propagator owning this stage
  float fFireFlushRatio   = 0;             ///< Ratio fired/flushed baskets to trigger basketizing
  int fId                 = 0;             ///< Unique stage id
  int fUserActionsStage   = 0;             ///< User actions stage to be executed right after
  int fFollowUpStage      = 0;             ///< In case there is a single follow-up store its id
  size_t fThrBasketCheck  = 0;             ///< Threshold for starting checking efficiency of basketizing
  size_t fNstaged         = 0;             ///< Total number of staged tracks
  size_t fNbasketized     = 0;             ///< Number of basketized handlers
  size_t fLocalHandlers   = 0;             ///< Local handlers (created only if basketizing is ON)
  bool fUniqueFollowUp    = false;         ///< All tracks go to single follow-up after this stage
  bool fEndStage          = false;         ///< Marker for stage at end of stepping
  bool fBasketized        = false;         ///< Stage is basketized
  Handlers_t fHandlers;                    ///< Array of handlers for the stage
  atomic_t<int> fCheckCountdown;           ///< Countdown fir checking basketizer efficiency

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  std::atomic_flag fCheckLock; ///< Lock for checking basketizers efficiency
#endif

private:
  VECCORE_ATT_HOST_DEVICE
  int CopyToFollowUps(Basket &output, TaskData *td);

  /** @brief  Check efficiency of basketizers. If less than threshold, flush and de-activate.
   * @return number of deactivated basketizers */
  VECCORE_ATT_HOST_DEVICE
  int CheckBasketizers(TaskData *td, size_t flush_threshold);

public:
  // The functions below are the interfaces for derived simulation stages.

  /** @brief Interface to create all handlers for the simulation stage
   *  @return Number of handlers created */
  VECCORE_ATT_HOST_DEVICE
  virtual int CreateHandlers() = 0;

  /** @brief Interface to select the handler matching a track */
  VECCORE_ATT_HOST_DEVICE
  virtual Handler *Select(Track *track, TaskData *td) = 0;

public:
  /** @brief Dummy SimulationStage constructor */
  VECCORE_ATT_HOST_DEVICE
  SimulationStage() {}

  /** @brief Standard SimulationStage constructor */
  VECCORE_ATT_HOST_DEVICE
  SimulationStage(ESimulationStage type, Propagator *prop);

  /** @brief Simulation stage destructor */
  VECCORE_ATT_HOST_DEVICE
  virtual ~SimulationStage();

  VECCORE_ATT_HOST_DEVICE
  SimulationStage(const SimulationStage &);

  VECCORE_ATT_HOST_DEVICE
  SimulationStage &operator=(const SimulationStage &);

  /** @brief Clone the stage and copy the existing handlers **/
  VECCORE_ATT_HOST_DEVICE
  virtual SimulationStage *Clone() const = 0;

  /** @brief Simulation stage name */
  VECCORE_ATT_HOST_DEVICE
  virtual const char *GetName() const = 0;

  /** @brief Get checking countdown */
  VECCORE_ATT_HOST_DEVICE
  int GetCheckCountdown() const
  {
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
    return fCheckCountdown;
#else
    return fCheckCountdown.load();
#endif
  }

  VECCORE_ATT_HOST_DEVICE
  void Info() const;

  //=== The stage processing methods === //

  /** @brief Set basketizing on/off */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetBasketizing(bool flag) { fBasketized = flag; }

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsBasketized() const { return fBasketized; }

  /** @brief Activate basketizing */
  VECCORE_ATT_HOST_DEVICE
  virtual void ActivateBasketizing(bool) {}

  /** @brief Count basketized handlers */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t CountBasketized()
  {
    for (int i = 0; i < GetNhandlers(); ++i)
      fNbasketized += size_t(fHandlers[i]->MayBasketize());
    return fNbasketized;
  }

  /** @brief Get number of basketized handlers */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t GetNbasketized() const { return fNbasketized; }

  /** @brief Count number of tracks stored in baskets */
  size_t GetNstored() const;

  /** @brief Process a basket of tracks marked for the stage
   *  @return Number of tracks processed
   */
  VECCORE_ATT_HOST_DEVICE
  int Process(TaskData *td);

  /** @brief Flush all tracks from the simulation stage basketizers and execute stage
   *  @return Number of tracks flushed
   */
  VECCORE_ATT_HOST_DEVICE
  int FlushAndProcess(TaskData *td);

  /** @brief Flush a handler and return the number of flushed tracks */
  VECCORE_ATT_HOST_DEVICE
  int FlushHandler(int i, TaskData *td, Basket &output);

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
  void SetFollowUpStage(ESimulationStage stage, bool unique = false)
  {
    fFollowUpStage  = (int)stage;
    fUniqueFollowUp = unique;
  }

  /** @brief Getter for type */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetFollowUpStage() const { return fFollowUpStage; }

  /** @brief Add next handler */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Handler *AddHandler(Handler *handler)
  {
    size_t id = fHandlers.size();
    handler->SetId(id);
    fHandlers.push_back(handler);
    fThrBasketCheck += handler->GetThreshold();
    fCheckCountdown = fThrBasketCheck;
    return handler;
  }

  /** @brief Getter for number of handlers */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNhandlers() const { return fHandlers.size(); }

  /** @brief Getter for a given handler */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Handler *GetHandler(int i) const { return fHandlers[i]; }

  /** @brief Getter for number of handlers */
  VECCORE_ATT_HOST_DEVICE
  bool HasLocalHandlers() const;

  /** @brief Getter for number of handlers */
  VECCORE_ATT_HOST_DEVICE
  size_t GetLocalHandlers() const { return fLocalHandlers; }

  /** @brief Getter for number of handlers */
  VECCORE_ATT_HOST_DEVICE
  void SetLocalHandlers(size_t mask) { fLocalHandlers = mask; }

  /** @brief Replace local handlers with LocalHandler */
  VECCORE_ATT_HOST_DEVICE
  virtual void ReplaceLocalHandlers();

  /** @brief Replace local handlers with LocalHandler */
  VECCORE_ATT_HOST_DEVICE
  void DeleteLocalHandlers();

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

  /** @brief Print handlers usage statistics */
  void PrintStatistics() const;
};

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant

#endif
