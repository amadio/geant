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

#if __cplusplus >= 201103L
#include <atomic>
#endif

class TTree;
class TFile;
class TStopwatch;
class TGeoVolume;
class PhysicsProcess;
class GeantTrack;
class GeantEvent;
class GeantBasket;
class GeantOutput;
class GeantBasketMgr;
class WorkloadManager;
class GeantThreadData;
class GeantVApplication;
class PrimaryGenerator;

class GeantPropagator : public TObject {
public:
  // data members to be made private
  Int_t fNthreads; /** Number of threads */
  Int_t fNevents;  /** Number of buffered */
  Int_t fNtotal;   /** Total number of events */
#if __cplusplus >= 201103L
  std::atomic<Long64_t> fNtransported; /** Number of transported tracks */
  std::atomic<Long64_t> fNprimaries;   /** Number of primary tracks */
  std::atomic<Long64_t> fNsafeSteps;   /** Number of fast steps within safety */
  std::atomic<Long64_t> fNsnextSteps;  /** Number of steps where full snext computation is needed */
  std::atomic<Long64_t> fNphysSteps;   /** Number of steps to physics process */
#endif
  Int_t fNprocesses;    /** Number of active processes */
  Int_t fElossInd;      /** Index of eloss process */
  Int_t fNstart;        /** Cumulated initial number of tracks */
  Int_t fMaxTracks;     /** Maximum number of tracks per event */
  Int_t fMaxThreads;    /** Maximum number of threads */
  Int_t fNminThreshold; /** Threshold for starting transporting a basket */
  Int_t fDebugTrk;      /** Track to debug */
  Int_t fMaxSteps;      /** Maximum number of steps per track */
  Int_t fNperBasket;    /** Number of tracks per basket */
  Int_t fMaxPerBasket;  /** Maximum number of tracks per basket */
  Int_t fMaxPerEvent;   /** Maximum number of tracks per event */
  Int_t fMaxDepth;      /** Maximum geometry depth */
  
  Double_t fMaxRes;     /** Maximum resident memory allowed [MBytes] */
  Double_t fNaverage;  /** Average number of tracks per event */
  Double_t fVertex[3]; /** Vertex position */
  Double_t fEmin;      /** Min energy threshold */
  Double_t fEmax;      /** Max energy threshold */
  Double_t fBmag;      /** Magnetic field */

  Bool_t fUsePhysics;       /** Enable/disable physics */
  Bool_t fUseDebug;         /** Use debug mode */
  Bool_t fUseGraphics;      /** Graphics mode */
  Bool_t fTransportOngoing; /** Flag for ongoing transport */
  Bool_t fSingleTrack;      /** Use single track transport mode */
  Bool_t fFillTree;         /** Enable I/O */
  Bool_t fUseMonitoring;    /** Monitoring thread */
  TMutex fTracksLock;       /** Mutex for adding tracks */

  WorkloadManager *fWMgr;          /** Workload manager */
  GeantVApplication *fApplication; /** User application */
  GeantOutput *fOutput;            /** Output object */
 
  TTree *fOutTree;    /** Output tree */
  TFile *fOutFile;    /** Output file */
  TStopwatch *fTimer; /** Timer */

  PhysicsProcess *fProcess; /** For now the only generic process pointing to the tabulated physics */
  PhysicsProcess *fVectorPhysicsProcess; /** interface to vector physics final state sampling */
  //   PhysicsProcess **fProcesses; //![fNprocesses] Array of processes
  GeantTrack_v *fStoredTracks; /** Stored array of tracks (history?) */
  PrimaryGenerator *fPrimaryGenerator;

  // Data per event
  Int_t *fNtracks;      /** ![fNevents] Number of tracks {array of [fNevents]} */
  GeantEvent **fEvents; /** ![fNevents]    Array of events */
  GeantThreadData **fThreadData; /** ![fNthreads] */

  static GeantPropagator *fgInstance;

public:

  /** @brief GeantPropagator constructor */
  GeantPropagator();

  /** @brief GeantPropagator destructor */
  virtual ~GeantPropagator();

#if __cplusplus >= 201103L

  /**
   * @brief Function that returned number of transported tracks (C++11)
   * @return Number of transported tracks
   */
  Long64_t GetNtransported() const { return fNtransported.load(); }
#else

  /**
   * @brief Function that returned number of transported tracks
   * @return Number of transported tracks
   */
  Long64_t GetNtransported() const { return fNtransported; }
#endif

  /**
   * @brief Function that returned temperary track
   * @details Temporary track for the current caller thread
   * 
   * @param tid Track ID 
   */
  GeantTrack &GetTempTrack(Int_t tid = -1);

  /**
   * @brief Function of addition track
   * 
   * @param track Track that should be added
   */
  Int_t AddTrack(GeantTrack &track);

  /**
   * @brief Function of dispatching track
   * 
   * @param track Track that should be dispatched
   */
  Int_t DispatchTrack(GeantTrack &track);

  /**
   * @brief  Function of stopping tracks
   * 
   * @param tracks Tracks that should be stopped
   * @param itr Track ID ?
   */
  void StopTrack(const GeantTrack_v &tracks, Int_t itr);

  /** @brief Function that returned index of elossity process */
  Int_t GetElossInd() const { return fElossInd; }

  /** @brief Function of loading geometry */
  Bool_t LoadGeometry(const char *filename = "geometry.root");
#if USE_VECGEOM_NAVIGATOR == 1

  /** @brief Function of loading VecGeom geometry */
  Bool_t LoadVecGeomGeometry();
#endif

  /**
   * @brief Function of importing tracks 
   * 
   * @param nevents Number of events
   * @param average Average value
   * @param startevent Start event
   * @param startslot Start slot
   */
  Int_t ImportTracks(Int_t nevents, Double_t average, Int_t startevent = 0, Int_t startslot = 0);
  
  /** @brief Initialization function */
  void Initialize();
  
  /**
   * @brief Instance function
   * 
   * @param ntotal Total number of tracks
   * @param nbuffered Number of buffered tracks
   */
  static GeantPropagator *Instance(Int_t ntotal = 0, Int_t nbuffered = 0);
  
  /**
   * @brief Propose step function
   * 
   * @param ntracks Number of threads
   * @param tracks Vector of tracks 
   * @param tid Track ID ?
   */
  void ProposeStep(Int_t ntracks, GeantTrack_v &tracks, Int_t tid);
  
  /**
   * @brief Apply Msc ????? (hard name)
   * 
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks
   * @param tid Track ID ?
   */
  void ApplyMsc(Int_t ntracks, GeantTrack_v &tracks, Int_t tid);
  //   PhysicsProcess  *Process(Int_t iproc) const {return fProcesses[iproc];}
  
  /**
   * @brief Function of definition of process
   * @return  Generic process pointing to the tabulated physics
   */
  PhysicsProcess *Process() const { return fProcess; }
  
  /**
   * @brief Function of geometry propagator 
   * 
   * @param geomfile Geometry file
   * @param nthreads Number of threads
   * @param graphics Graphics (by default False)
   * @param single ????? (by default False)
   */
  void PropagatorGeom(const char *geomfile = "geometry.root", Int_t nthreads = 4,
                      Bool_t graphics = kFALSE, Bool_t single = kFALSE);

private:

  /** @brief Not implemented */
  GeantPropagator(const GeantPropagator &);
  
  /** @brief Not implemented */
  GeantPropagator &operator=(const GeantPropagator &);

  ClassDef(GeantPropagator, 1)
};
#endif
