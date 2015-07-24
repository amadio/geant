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

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef GEANT_TRACK
#include "GeantTrack.h"
#endif

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

#include <vector>
#include <atomic>

#include "Geant/Typedefs.h"

class TTree;
class TFile;
class TStopwatch;
class PhysicsProcess;
class GeantEvent;
class GeantBasket;
class GeantBasketMgr;
class WorkloadManager;
class GeantVApplication;
class PrimaryGenerator;
class TaskBroker;

#include "GeantFwd.h"

class GeantPropagator : public TObject {
public:
  /**
   * @brief Monitoring type
   */
  enum EGeantMonitoringType {
    kMonQueue = 0,
    kMonMemory,
    kMonBasketsPerVol,
    kMonVectors,
    kMonConcurrency,
    kMonTracksPerEvent,
    kMonTracks
  };
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;
  // data members to be made private
  Int_t fNthreads;                     /** Number of worker threads */
  Int_t fNevents;                      /** Number of buffered events */
  Int_t fNtotal;                       /** Total number of events to be transported */
  std::atomic<Long64_t> fNtransported; /** Number of transported tracks */
  std::atomic<Long64_t> fNprimaries;   /** Number of primary tracks */
  std::atomic<Long64_t> fNsafeSteps;   /** Number of fast steps within safety */
  std::atomic<Long64_t> fNsnextSteps;  /** Number of steps where full snext computation is needed */
  std::atomic<Long64_t> fNphysSteps;   /** Number of steps to physics process */
  std::atomic_flag fFeederLock;        /** Atomic flag to protect the particle feeder */
  std::atomic_int fPriorityEvents;     /** Number of prioritized events */
  BitSet *fDoneEvents;                 /** Array of bits marking done events */
  Int_t fNprocesses;                   /** Number of active physics processes */
  Int_t fNstart;                       /** Cumulated initial number of tracks */
  Int_t fMaxTracks;                    /** Maximum number of tracks per event */
  Int_t fMaxThreads;                   /** Maximum number of threads */
  Int_t fNminThreshold;                /** Threshold for starting transporting a basket */
  Int_t fDebugEvt;                     /** Event to debug */
  Int_t fDebugTrk;                     /** Track to debug */
  Int_t fDebugStp;                     /** Step to start debugging */
  Int_t fDebugRep;                     /** Number of steps to debug */
  Int_t fMaxSteps;                     /** Maximum number of steps per track */
  Int_t fNperBasket;                   /** Number of tracks per basket */
  Int_t fMaxPerBasket;                 /** Maximum number of tracks per basket */
  Int_t fMaxPerEvent;                  /** Maximum number of tracks per event */
  Int_t fMaxDepth;                     /** Maximum geometry depth */
  Int_t fLearnSteps;                   /** Number of steps needed for the learning phase */
  Int_t fLastEvent;                    /** Last transported event */
  Float_t fPriorityThr;                /** Threshold for prioritizing events */

  Double_t fMaxRes;    /** Maximum resident memory allowed [MBytes] */
  Double_t fNaverage;  /** Average number of tracks per event */
  Double_t fVertex[3]; /** Vertex position */
  Double_t fEmin;      /** Min energy threshold */
  Double_t fEmax;      /** Max energy threshold */
  Double_t fBmag;      /** Magnetic field */

  Bool_t fUsePhysics;       /** Enable/disable physics */
  Bool_t fUseDebug;         /** Use debug mode */
  Bool_t fUseGraphics;      /** Graphics mode */
  Bool_t fUseStdScoring;    /** Use standard scoring */
  Bool_t fTransportOngoing; /** Flag for ongoing transport */
  Bool_t fSingleTrack;      /** Use single track transport mode */
  Bool_t fUseMonitoring;    /** Monitoring different features */
  Bool_t fUseAppMonitoring; /** Monitoring the application */
  TMutex fTracksLock;       /** Mutex for adding tracks */

  WorkloadManager *fWMgr;             /** Workload manager */
  GeantVApplication *fApplication;    /** User application */
  GeantVApplication *fStdApplication; /** Standard application */

  TStopwatch *fTimer; /** Timer */

  PhysicsProcess *fProcess;              /** For now the only generic process pointing to the tabulated physics */
  PhysicsProcess *fVectorPhysicsProcess; /** interface to vector physics final state sampling */
  //   PhysicsProcess **fProcesses; //![fNprocesses] Array of processes
  GeantTrack_v *fStoredTracks;         /** Stored array of tracks (history?) */
  PrimaryGenerator *fPrimaryGenerator; /** Primary generator */

  // Data per event
  Int_t *fNtracks;             /** ![fNevents] Number of tracks {array of [fNevents]} */
  GeantEvent **fEvents;        /** ![fNevents]    Array of events */
  GeantTaskData **fThreadData; /** ![fNthreads] Data private to threads */

  static GeantPropagator *fgInstance;

  /** @brief Initialization function */
  void Initialize();

  /** @brief Initialization function */
  void InitializeAfterGeom();

  /** @brief Function for loading geometry */
  Bool_t LoadGeometry(const char *filename = "geometry.root");
#if USE_VECGEOM_NAVIGATOR == 1

  /** @brief Function for loading VecGeom geometry */
  Bool_t LoadVecGeomGeometry();
#endif

public:
  /** @brief GeantPropagator constructor */
  GeantPropagator();

  /** @brief GeantPropagator destructor */
  virtual ~GeantPropagator();

  /**
   * @brief Function that returns the number of prioritized events (C++11)
   * @return Number of prioritized events
   */
  Int_t GetNpriority() const { return fPriorityEvents.load(); }

  /**
   * @brief Function that returns the number of transported tracks (C++11)
   * @return Number of transported tracks
   */
  Long64_t GetNtransported() const { return fNtransported.load(); }

  /**
   * @brief Function that returns a temporary track object per thread
   * @details Temporary track for the current caller thread
   *
   * @param tid Track ID
   */
  GeantTrack &GetTempTrack(Int_t tid = -1);

  /**
   * @brief Function to add a track to the scheduler
   *
   * @param track Track that should be added
   */
  Int_t AddTrack(GeantTrack &track);

  /**
   * @brief Function to dispatch a track
   *
   * @param track Track that should be dispatched
   */
  Int_t DispatchTrack(GeantTrack &track, GeantTaskData *td);

  /**
   * @brief  Function for marking a track as stopped
   *
   * @param tracks Track array container
   * @param itr Track id
   */
  void StopTrack(const GeantTrack_v &tracks, Int_t itr);

  /**
   * @brief Feeder for importing tracks
   *
   * @param td Thread data object
   * @param init Flag specifying if this is the first call
   * @return Number of injected baskets
   */
  Int_t Feeder(GeantTaskData *td);

  /** @brief Check if transport is feeding with new tracks. */
  inline Bool_t IsFeeding() {
    Bool_t feeding = fFeederLock.test_and_set(std::memory_order_acquire);
    if (feeding)
      return kTRUE;
    fFeederLock.clear(std::memory_order_release);
    return kFALSE;
  }

  /**
   * @brief Function for importing tracks
   *
   * @param nevents Number of events
   * @param average Average number of tracks
   * @param startevent Start event
   * @param startslot Start slot
   */
  Int_t ImportTracks(Int_t nevents, Int_t startevent, Int_t startslot, GeantTaskData *td);

  /**
   * @brief Instance function returning the singleton pointer
   *
   * @param ntotal Total number of tracks
   * @param nbuffered Number of buffered tracks
   */
  static GeantPropagator *Instance(Int_t ntotal = 0, Int_t nbuffered = 0, Int_t nthreads = 0);

  /**
   * @brief Propose the physics step for an array of tracks
   *
   * @param ntracks Number of ttracks
   * @param tracks Vector of tracks
   * @param td Thread data
   */
  void ProposeStep(Int_t ntracks, GeantTrack_v &tracks, GeantTaskData *td);

  /**
   * @brief Apply multiple scattering process
   *
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks
   * @param td Thread data
   */
  void ApplyMsc(Int_t ntracks, GeantTrack_v &tracks, GeantTaskData *td);
  //   PhysicsProcess  *Process(Int_t iproc) const {return fProcesses[iproc];}

  /**
   * @brief Getter for the process
   * @return  Generic process pointing to the tabulated physics
   */
  PhysicsProcess *Process() const { return fProcess; }

  /**
   * @brief Entry point to start simulation with GeantV
   *
   * @param geomfile Geometry file
   * @param nthreads Number of threads
   * @param graphics Graphics (by default False)
   * @param single Transport single tracks rather than vectors (by default False)
   */
  void PropagatorGeom(const char *geomfile = "geometry.root", Int_t nthreads = 4, Bool_t graphics = kFALSE,
                      Bool_t single = kFALSE);

  /** @brief Function returning the number of monitored features */
  int GetMonFeatures() const;

  /** @brief Check if a monitoring feature is enabled */
  bool IsMonitored(GeantPropagator::EGeantMonitoringType feature) const;

  /** @brief Enable monitoring a feature */
  void SetMonitored(EGeantMonitoringType feature, bool flag = true);

  /** @brief Setter for the global transport threshold */
  void SetNminThreshold(int thr);

  /** @brief  Getter for task broker */
  TaskBroker *GetTaskBroker();

  /** @brief  Setter for task broker */
  void SetTaskBroker(TaskBroker *broker);

  /** @brief Function checking if transport is completed */
  Bool_t TransportCompleted() const { return ((Int_t)fDoneEvents->FirstNullBit() >= fNtotal); }

private:
  /** @brief Copy constructor not implemented */
  GeantPropagator(const GeantPropagator &);

  /** @brief Assignment operator not implemented */
  GeantPropagator &operator=(const GeantPropagator &);

  ClassDef(GeantPropagator, 1)
};
#endif
