//===--- Propagator.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Propagator.h
 * @brief Implementation of propogator in Geant-V prototype.
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PROPAGATOR
#define GEANT_PROPAGATOR

#include <vector>
#include <atomic>
#include <mutex>

#include "Geant/Typedefs.h"
#ifdef USE_ROOT
class TTree;
class TFile;
class TStopwatch;
#else
#include "base/Stopwatch.h"
#endif

#include "GeantConfig.h"
#include "Geant/Track.h"

class PhysicsInterface;

class GUFieldPropagator;
class VVectorField;

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class RunManager;
class Event;
class UserApplication;
class GeantBasket;
class GeantBasketMgr;
class WorkloadManager;
class PrimaryGenerator;
class MCTruthMgr;
class TaskBroker;
class SimulationStage;
class TrackManager;
// class VVectorField;
// class GUFieldPropagator;

// #include "Geant/Fwd.h"

class Propagator {
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  // On cuda there is one propagator per thread.  So (for now), no need
  // for atomics.
  template <typename T>
  using atomic_t = T;
#else
  template <typename T>
  using atomic_t             = std::atomic<T>;
#endif

public:
  GeantConfig *fConfig = nullptr; /** Run configuration*/
  RunManager *fRunMgr  = nullptr; /** Run manager */
  int fNuma            = -1;      /** NUMA id */

  int fNthreads = 0;            /** Number of worker threads */
  int fNbuff    = 0;            /** Number of buffered events */
  int fNtotal   = 0;            /** Total number of events to be transported */
  atomic_t<long> fNtransported; /** Number of transported tracks */
  atomic_t<long> fNsteps;       /** Total number of steps */
  atomic_t<long> fNsnext;       /** Total number of calls to getting distance to next boundary */
  atomic_t<long> fNphys;        /** Total number of steps to physics processes */
  atomic_t<long> fNmag;         /** Total number of partial steps in magnetic field */
  atomic_t<long> fNsmall;       /** Total number of small steps taken */
  atomic_t<long> fNcross;       /** Total number of boundaries crossed */
  atomic_t<long> fNpushed;      /** Total number of pushes of 1E-3 */
  atomic_t<long> fNkilled;      /** Total number of killed tracks */
  atomic_t<int> fNbfeed;        /** Number of baskets fed from server */

  bool fTransportOngoing = false; /** Flag for ongoing transport */
  bool fSingleTrack      = false; /** Use single track transport mode */
  bool fDoFastSim        = false; /** Include Fast Sim stage */

  WorkloadManager *fWMgr           = nullptr; /** Workload manager */
  UserApplication *fApplication    = nullptr; /** User application */
  UserApplication *fStdApplication = nullptr; /** Standard application */

#ifdef USE_ROOT
  TStopwatch *fTimer = nullptr; /** Timer */
#else
  vecgeom::Stopwatch *fTimer = nullptr; /** Timer */
#endif

  PhysicsInterface *fPhysicsInterface = nullptr; /** The new, real physics interface */
  PrimaryGenerator *fPrimaryGenerator = nullptr; /** Primary generator */
  MCTruthMgr *fTruthMgr               = nullptr; /** MCTruth manager */
  TrackManager *fTrackMgr             = nullptr; /** Track manager */
  vector_t<SimulationStage *> fStages;           /** Simulation stages */

  // Data per event
  int *fNtracks     = nullptr; /** ![fNbuff] Number of tracks per slot */
  Event **fEvents   = nullptr; /** ![fNbuff]    Array of events */
  bool fCompleted   = false;   /** Completion flag */
  bool fInitialFeed = false;   /** Flag marking that events were injected */
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  std::mutex fStopperLock; /** Mutex for the stop operation */
#endif

  /** @brief Initialization function */
  void Initialize();

  /** @brief Initialization function */
  void InitializeAfterGeom();

public:
  /** @brief Propagator constructor
   * @param ntotal Total number of tracks
   * old param nbuffered Number of buffered tracks
   */
  VECCORE_ATT_HOST_DEVICE
  Propagator(int nthreads = 0);

  /** @brief Copy constructor */
  Propagator(const Propagator &);

  /** @brief Propagator destructor */
  virtual ~Propagator();

  /**
   * @brief Function that returns the number of transported tracks (C++11)
   * @return Number of transported tracks
   */
  GEANT_FORCE_INLINE
  long GetNtransported() const { return fNtransported; }

  /** @brief Stop the transport threads */
  void StopTransport();

  /**
   * @brief  Function for marking a track as stopped
   *
   * @param track Track to be stopped
   * @param itr Track id
   */
  void StopTrack(Track *track, TaskData *td);

  /**
   * @brief Setter for the real physics interface
   * @param[in] physics Pointer to the real physics framework object.
   */
  void SetPhysicsInterface(PhysicsInterface *physics) { fPhysicsInterface = physics; }

  /**
   * @brief   Getter for the real physics interface
   * @return  Generic physics interface pointing to the real physics framework
   */
  PhysicsInterface *GetPhysicsInterface() const { return fPhysicsInterface; }

  /** @brief Entry point to start simulation with GeantV */
  static void RunSimulation(Propagator *prop, int nthreads);

  /**
   * @brief Entry point to start simulation with GeantV
   *
   * @param geomfile Geometry file
   * @param nthreads Number of threads
   * @param graphics Graphics (by default False)
   * @param single Transport single tracks rather than vectors (by default False)
   * @param Execution using TBB tasks instead of static threads (by default False)
   */
  void PropagatorGeom(int nthreads);

  ///** @brief Function returning the number of monitored features */
  // int GetMonFeatures() const;

  ///** @brief Setter for the global transport threshold */
  // void SetNminThreshold(int thr);

  /** @brief  Getter for task broker */
  TaskBroker *GetTaskBroker();

  /** @brief  Setter for task broker */
  void SetTaskBroker(TaskBroker *broker);

  /** @brief  Synchronize with run configuration */
  void SetConfig(GeantConfig *config);

  /** @brief  Register a simulation stage */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int RegisterStage(SimulationStage *stage)
  {
    fStages.push_back(stage);
    return (fStages.size() - 1);
  }

  /** @brief  Inspect simulation stages */
  VECCORE_ATT_HOST_DEVICE
  void InspectStages() const;

  /** @brief  Getter for a simulation stage */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  SimulationStage *GetStage(ESimulationStage id) const { return fStages[int(id)]; }

  /** @brief  Getter for the number of simulation stages */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNstages() const { return fStages.size(); }

  /** @brief Function creating all simulation stages for a propagator */
  int CreateSimulationStages();

  /** @brief Function allowing to retrieve the next simulation stage for a track */
  VECCORE_ATT_HOST_DEVICE
  int GetNextStage(Track &track, int current) const;

  /** @brief Setter for locality */
  void SetNuma(int numa);

private:
  /** @brief Assignment operator not implemented */
  Propagator &operator=(const Propagator &);
};

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant

#endif
