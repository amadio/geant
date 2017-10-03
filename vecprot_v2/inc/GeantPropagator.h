//===--- GeantPropagator.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantPropagator.h
 * @brief Implementation of propogator in Geant-V prototype.
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PROPAGATOR
#define GEANT_PROPAGATOR

#ifndef GEANT_TRACK
#include "GeantTrackVec.h"
#endif

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
#include "base/BitSet.h"
using veccore::BitSet;

#include "GeantConfig.h"

class PhysicsInterface;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTrack_v;
class GeantRunManager;
class GeantEvent;
class GeantVApplication;
class PhysicsProcessOld;
class GeantBasket;
class GeantBasketMgr;
class WorkloadManager;
class GeantVTaskMgr;
class PrimaryGenerator;
class MCTruthMgr;
class TaskBroker;
class SimulationStage;
class TrackManager;

class GeantPropagator {
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  // On cuda there is one propagator per thread.  So (for now), no need
  // for atomics.
  template <typename T>
  using atomic_t = T;
#else
  template <typename T>
  using atomic_t = std::atomic<T>;
#endif

public:
  GeantConfig *fConfig = nullptr;      /** Run configuration*/
  GeantRunManager *fRunMgr = nullptr;  /** Run manager */
  int fNuma = -1;                      /** NUMA id */

  int fNthreads = 0;                   /** Number of worker threads */
  int fNbuff = 0;                      /** Number of buffered events */
  int fNtotal = 0;                     /** Total number of events to be transported */
  atomic_t<long> fNtransported;        /** Number of transported tracks */
  atomic_t<long> fNsteps;              /** Total number of steps */
  atomic_t<long> fNsnext;              /** Total number of calls to getting distance to next boundary */
  atomic_t<long> fNphys;               /** Total number of steps to physics processes */
  atomic_t<long> fNmag;                /** Total number of partial steps in magnetic field */
  atomic_t<long> fNsmall;              /** Total number of small steps taken */
  atomic_t<long> fNcross;              /** Total number of boundaries crossed */
  atomic_t<long> fNpushed;             /** Total number of pushes of 1E-3 */
  atomic_t<long> fNkilled;             /** Total number of killed tracks */
  atomic_t<int>  fNidle;               /** Number of idle threads */
  atomic_t<int>  fNbfeed;              /** Number of baskets fed from server */

  bool fTransportOngoing = false;      /** Flag for ongoing transport */
  bool fSingleTrack = false;           /** Use single track transport mode */

  WorkloadManager   *fWMgr = nullptr;           /** Workload manager */
  GeantVApplication *fApplication = nullptr;    /** User application */
  GeantVApplication *fStdApplication = nullptr; /** Standard application */
  GeantVTaskMgr     *fTaskMgr = nullptr;        /** GeantV task manager */

  #ifdef USE_ROOT
  TStopwatch *fTimer = nullptr;                 /** Timer */
  #else
  vecgeom::Stopwatch *fTimer = nullptr;         /** Timer */
  #endif

  PhysicsProcessOld *fProcess = nullptr;              /** For now the only generic process pointing to the tabulated physics */
  PhysicsProcessOld *fVectorPhysicsProcess = nullptr; /** interface to vector physics final state sampling */
  PhysicsInterface *fPhysicsInterface = nullptr;     /** The new, real physics interface */
  GeantTrack_v *fStoredTracks = nullptr;         /** Stored array of tracks (history?) */
  PrimaryGenerator *fPrimaryGenerator = nullptr; /** Primary generator */
  MCTruthMgr *fTruthMgr = nullptr;               /** MCTruth manager */
  TrackManager *fTrackMgr = nullptr;             /** Track manager */
  vector_t<SimulationStage *> fStages;           /** Simulation stages */

  // Data per event
  int *fNtracks = nullptr;        /** ![fNbuff] Number of tracks per slot */
  GeantEvent **fEvents = nullptr; /** ![fNbuff]    Array of events */
  bool fCompleted = false;     /** Completion flag */
  bool fInitialFeed = false;   /** Flag marking that events were injected */
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  std::mutex fStopperLock;     /** Mutex for the stop operation */
#endif

  /** @brief Initialization function */
  void Initialize();

  /** @brief Initialization function */
  void InitializeAfterGeom();

public:
  /** @brief GeantPropagator constructor
   * @param ntotal Total number of tracks
   * old param nbuffered Number of buffered tracks
   */
  VECCORE_ATT_HOST_DEVICE
  GeantPropagator(int nthreads = 0);

  /** @brief Copy constructor */
  GeantPropagator(const GeantPropagator &);

  /** @brief GeantPropagator destructor */
  virtual ~GeantPropagator();

  /**
   * @brief Function that returns the number of transported tracks (C++11)
   * @return Number of transported tracks
   */
  GEANT_FORCE_INLINE
  long GetNtransported() const { return fNtransported; }

  /** @brief Check if the propagator threads are frozen */
  bool IsIdle() const;

  /** @brief Get number of idle workers */
  GEANT_FORCE_INLINE
  int GetNidle() const { return fNidle; }

  /** @brief Get number of pending baskets */
  GEANT_FORCE_INLINE
  int GetNpending() const;

  /** @brief Get number of working threads */
  GEANT_FORCE_INLINE
  int GetNworking() const { return (fNthreads - fNidle); }

  /** @brief Stop the transport threads */
  void StopTransport();

  /**
   * @brief Function to add a track to the scheduler
   *
   * @param track Track that should be added
   */
  int AddTrack(GeantTrack &track);

  /**
   * @brief Function to dispatch a track
   *
   * @param track Track that should be dispatched
   */
  int DispatchTrack(GeantTrack &track, GeantTaskData *td);

  /**
   * @brief  Function for marking a track as stopped
   *
   * @param track Track to be stopped
   * @param itr Track id
   */
  void StopTrack(GeantTrack *track, GeantTaskData *td);

  /**
   * @brief  Function for marking a track as stopped
   *
   * @param tracks Track array container
   * @param itr Track id
   */
  void StopTrack(const GeantTrack_v &tracks, int itr, GeantTaskData *td);

  /**
   * @brief Propose the physics step for an array of tracks
   *
   * @param ntracks Number of ttracks
   * @param tracks Vector of tracks
   * @param td Thread data
   */
  void ProposeStep(int ntracks, GeantTrack_v &tracks, GeantTaskData *td);

  /**
   * @brief Apply multiple scattering process
   *
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks
   * @param td Thread data
   */
  void ApplyMsc(int ntracks, GeantTrack_v &tracks, GeantTaskData *td);
  //   PhysicsProcess  *Process(int iproc) const {return fProcesses[iproc];}

  /**
   * @brief Getter for the process
   * @return  Generic process pointing to the tabulated physics
   */
  VECCORE_ATT_HOST_DEVICE
  PhysicsProcessOld *Process() const { return fProcess; }

  /**
   * @brief Setter for the real physics interface
   * @param[in] physics Pointer to the real physics framework object.
   */
  void SetPhysicsInterface(PhysicsInterface* physics) { fPhysicsInterface = physics; }

  /**
   * @brief   Getter for the real physics interface
   * @return  Generic physics interface pointing to the real physics framework
   */
  PhysicsInterface *GetPhysicsInterface() const { return fPhysicsInterface; }


  /** @brief Entry point to start simulation with GeantV */
  static void RunSimulation(GeantPropagator *prop, int nthreads);

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
  //int GetMonFeatures() const;

  ///** @brief Setter for the global transport threshold */
  //void SetNminThreshold(int thr);

  /** @brief  Getter for task broker */
  TaskBroker *GetTaskBroker();

  /** @brief  Setter for task broker */
  void SetTaskBroker(TaskBroker *broker);

  /** @brief  Synchronize with run configuration */
  void SetConfig(GeantConfig* config);

  /** @brief  Share work with some other propagator */
  int ShareWork(GeantPropagator &other);

  /** @brief  Register a simulation stage */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int RegisterStage(SimulationStage *stage)
  {
    fStages.push_back(stage);
    return ( fStages.size() - 1);
  }

  /** @brief  Inspect simulation stages */
  VECCORE_ATT_HOST_DEVICE
  void InspectStages();

  /** @brief  Getter for a simulation stage */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  SimulationStage *GetStage(ESimulationStage id) { return fStages[int(id)]; }

  /** @brief  Getter for the number of simulation stages */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNstages() { return fStages.size(); }

  /** @brief Function creating all simulation stages for a propagator */
  #ifdef USE_REAL_PHYSICS
    int CreateSimulationStages();
  #else
    VECCORE_ATT_HOST_DEVICE
    int CreateSimulationStages();
  #endif

  /** @brief Function allowing to retrieve the next simulation stage for a track */
  VECCORE_ATT_HOST_DEVICE
  int GetNextStage(GeantTrack &track, int current);

  /** @brief Setter for locality */
  void SetNuma(int numa);

private:
  /** @brief Assignment operator not implemented */
  GeantPropagator &operator=(const GeantPropagator &);

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
