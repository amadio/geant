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
class PhysicsProcessOld;
class PhysicsInterface;
class GeantRunManager;
class GeantEvent;
class GeantBasket;
class GeantBasketMgr;
class WorkloadManager;
class GeantVApplication;
class GeantVTaskMgr;
class PrimaryGenerator;
class MCTruthMgr;
class TaskBroker;

#include "GeantFwd.h"
#include "GeantConfig.h"

class GeantPropagator {

public:
  GeantConfig *fConfig;                            /** Run configuration*/
  GeantRunManager *fRunMgr;                        /** Run manager */

  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;
  // data members to be made private
  int fNthreads;                                   /** Number of worker threads */
  int fNbuff;                                      /** Number of buffered events */
  int fNtotal;                                     /** Total number of events to be transported */
  std::atomic<long> fNtransported;                 /** Number of transported tracks */
  std::atomic<long> fNprimaries;                   /** Number of primary tracks */
  std::atomic<long> fNsteps;                       /** Total number of steps */
  std::atomic<long> fNsnext;                       /** Total number of calls to getting distance to next boundary */
  std::atomic<long> fNphys;                        /** Total number of steps to physics processes */
  std::atomic<long> fNmag;                         /** Total number of partial steps in magnetic field */
  std::atomic<long> fNsmall;                       /** Total number of small steps taken */
  std::atomic<long> fNcross;                       /** Total number of boundaries crossed */

  bool fTransportOngoing;      /** Flag for ongoing transport */
  bool fSingleTrack;           /** Use single track transport mode */
    
  WorkloadManager *fWMgr;             /** Workload manager */
  GeantVApplication *fApplication;    /** User application */
  GeantVApplication *fStdApplication; /** Standard application */
  GeantVTaskMgr     *fTaskMgr;        /** GeantV task manager */

  #ifdef USE_ROOT
  TStopwatch *fTimer; /** Timer */
  #else
  vecgeom::Stopwatch *fTimer; /** Timer */
  #endif

  PhysicsProcessOld   *fProcess;              /** For now the only generic process pointing to the tabulated physics */
  PhysicsProcessOld   *fVectorPhysicsProcess; /** interface to vector physics final state sampling */
  PhysicsInterface *fPhysicsInterface;     /** The new, real physics interface */

  //   PhysicsProcess **fProcesses; //![fNprocesses] Array of processes
  GeantTrack_v *fStoredTracks;         /** Stored array of tracks (history?) */
  PrimaryGenerator *fPrimaryGenerator; /** Primary generator */
  MCTruthMgr *fTruthMgr;               /** MCTruth manager */

  // Data per event
  int *fNtracks;               /** ![fNbuff] Number of tracks per slot */
  GeantEvent **fEvents;        /** ![fNbuff]    Array of events */
  GeantTaskData **fThreadData; /** ![fNthreads] Data private to threads */

  /** @brief Initialization function */
  void Initialize();

  /** @brief Initialization function */
  void InitializeAfterGeom();
  
public:
  /** @brief GeantPropagator constructor */
  GeantPropagator();

  /** @brief GeantPropagator destructor */
  virtual ~GeantPropagator();

  /**
   * @brief Function that returns the number of transported tracks (C++11)
   * @return Number of transported tracks
   */
  long GetNtransported() const { return fNtransported.load(); }

  /**
   * @brief Function that returns the number of primary tracks
   * @return Number of primary tracks
   */
  long GetNprimaries() const { return fNprimaries.load(); }

  /**
   * @brief Function that returns a temporary track object per thread
   * @details Temporary track for the current caller thread
   *
   * @param tid Track ID
   */
  GeantTrack &GetTempTrack(int tid = -1);

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
   * @param tracks Track array container
   * @param itr Track id
   */
  void StopTrack(const GeantTrack_v &tracks, int itr);

  /**
   * @brief Instance function returning the singleton pointer
   *
   * @param ntotal Total number of tracks
   * @param nbuffered Number of buffered tracks
   */
  VECCORE_ATT_HOST_DEVICE
  static GeantPropagator *NewInstance(int nthreads = 0);

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


  /**
   * @brief Entry point to start simulation with GeantV
   *
   * @param geomfile Geometry file
   * @param nthreads Number of threads
   * @param graphics Graphics (by default False)
   * @param single Transport single tracks rather than vectors (by default False)
   * @param Execution using TBB tasks instead of static threads (by default False)
   */
  void PropagatorGeom(const char *geomfile = "geometry.root", int nthreads = 4, bool graphics = false,
                      bool single = false);

  ///** @brief Function returning the number of monitored features */
  //int GetMonFeatures() const;

  ///** @brief Setter for the global transport threshold */
  //void SetNminThreshold(int thr);

  /** @brief  Getter for task broker */
  TaskBroker *GetTaskBroker();

  /** @brief  Setter for task broker */
  void SetTaskBroker(TaskBroker *broker);

  void SetConfig(GeantConfig* config);

private:
  /** @brief Copy constructor not implemented */
  GeantPropagator(const GeantPropagator &);

  /** @brief Assignment operator not implemented */
  GeantPropagator &operator=(const GeantPropagator &);

};
#endif
