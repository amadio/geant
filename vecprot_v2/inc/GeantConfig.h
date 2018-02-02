#ifndef GEANT_RUN_CONFIG_H
#define GEANT_RUN_CONFIG_H

#ifndef VECCORE_CUDA
#include <string>
#endif

class GeantVApplication;
class GeantVTaskMgr;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantConfig
{
public:
	/**
	* @brief Run mode
	*/
	enum ERunMode {
    kGenerator = 0,
    kExternalLoop,
    kHPC
  };

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

  ERunMode fRunMode = kGenerator; /** GeantV run mode */
  int fNtotal = 1000;     /** Total number of events to be transported */
  int fNbuff = 100;       /** Number of buffered events */
  int fNprocesses = 3;    /** Number of active physics processes */
  int fNstart = 0;        /** Cumulated initial number of tracks */
  int fMaxTracks = 0;     /** Maximum number of tracks per event */
  int fMaxThreads = 1000; /** Maximum number of threads */
  int fNminThreshold = 16;/** Threshold for starting transporting a basket */
  int fNvecThreshold = 8; /** Threshold for executing in vector mode */
  int fDebugEvt = -1;     /** Event to debug */
  int fDebugTrk = -1;     /** Track to debug */
  int fDebugStp = -1;     /** Step to start debugging */
  int fDebugRep = -1;     /** Number of steps to debug */
  int fMaxSteps = 10000;  /** Maximum number of steps per track */
  int fNperBasket = 16;   /** Number of tracks per basket */
  int fMaxPerBasket = 256;/** Maximum number of tracks per basket */
  int fMaxPerEvent = 0;   /** Maximum number of tracks per event */
  int fMaxDepth = 0;      /** Maximum geometry depth */
  int fLearnSteps = 0;    /** Number of steps needed for the learning phase */
  int fLastEvent = 0;     /** Last transported event */
  float fPriorityThr = 0; /** Threshold for prioritizing events */
  int fNstepsKillThr = 50000; /** Threshold in number of steps to kill a track */
  int fNminReuse = 10000; /** Minimum number of transported tracks to be reused without re-basketizing */
  int fNstackLanes = 10;  /** Number of stacked lanes in the stack-like buffers */
  int fNmaxBuffSpill = 128; /** Maximum number of tracks spilled from the stack-like buffer per stepping iteration */

  double fMaxRes = 0;     /** Maximum resident memory allowed [MBytes] */
  double fMaxVirt = 0;    /** Maximum virtual memory allowed [MBytes] */
  double fNaverage = 0;   /** Average number of tracks per event */
  double fVertex[3];      /** Vertex position */
  double fEmin = 1.E-4;   /** Min energy threshold [GeV] */
  double fEmax = 10;      /** Max energy threshold [GeV] */
  // double fBfieldMag = 0.0;   /** Magnitude of field in case of const field [kiloGauss] */
  double fEpsilonRK = 3.e-4;  /** Relative error in RK integration */
  float fFireFlushRatio = 2.; /** Ratio fired/flush to trigger basketizing */

  bool fUsePhysics = true;   /** Enable/disable physics */
  bool fUseRungeKutta = false; /** Enable/disable Runge-Kutta integration in field */
  bool fUseDebug = false; /** Use debug mode */
  bool fUseGraphics = false;   /** Graphics mode */
  bool fUseStdScoring = false; /** Use standard scoring */
  bool fUseV3 = false;         /** Use version 3 of the scheduler */
  bool fUseNuma = false;       /** Use NUMA */
  bool fUseVectorizedGeom = false; /** Use vectorized geometry */

  int fMonQueue = 0;      /** Monitor the work queue */
  int fMonMemory = 0;     /** Monitor the memory */
  int fMonBasketsPerVol = 0; /** Monitor baskets per volume */
  int fMonVectors = 0;    /** Monitor vector scheduling */
  int fMonConcurrency = 0;/** Monitor concurrency */
  int fMonTracksPerEvent = 0; /** Monitor tracks status per event */
  int fMonTracks = 0;     /** Monitor number of tracks */

  bool fFillTree = false; /** Enable I/O */
  bool fUseMonitoring = false; /** Monitoring different features */
  bool fUseAppMonitoring = false; /** Monitoring the application */
  int  fTreeSizeWriteThreshold = 100000; /** Maximum size of the tree (before automatic writing) **/
  bool fConcurrentWrite = true;  /** Switch between single and mutlithreaded writing */
#ifndef VECCORE_CUDA
  std::string fGeomFileName; /** Geometry file name */
#endif

public:
  VECCORE_ATT_DEVICE
  GeantConfig() {};

  VECCORE_ATT_DEVICE
  ~GeantConfig() {}

  /** @brief Check if a monitoring feature is enabled */
  bool IsMonitored(GeantConfig::EGeantMonitoringType feature) const {
  // Check if a given feature is monitored
    switch (feature) {
      case GeantConfig::kMonQueue:
      return fMonQueue;
      case GeantConfig::kMonMemory:
        return fMonMemory;
      case GeantConfig::kMonBasketsPerVol:
        return fMonBasketsPerVol;
      case GeantConfig::kMonVectors:
        return fMonVectors;
      case GeantConfig::kMonConcurrency:
        return fMonConcurrency;
      case GeantConfig::kMonTracksPerEvent:
        return fMonTracksPerEvent;
      case GeantConfig::kMonTracks:
        return fMonTracks;
    }
    return false;
  }

  /** @brief Enable monitoring a feature */
  void SetMonitored(GeantConfig::EGeantMonitoringType feature, bool flag = true) {
  // Enable/disable monitoring for a feature
    int value = (int)flag;
    switch (feature) {
    case GeantConfig::kMonQueue:
      fMonQueue = value;
      break;
    case GeantConfig::kMonMemory:
      fMonMemory = value;
      break;
    case GeantConfig::kMonBasketsPerVol:
      fMonBasketsPerVol = value;
      break;
    case GeantConfig::kMonVectors:
      fMonVectors = value;
      break;
    case GeantConfig::kMonConcurrency:
      fMonConcurrency = value;
      break;
    case GeantConfig::kMonTracksPerEvent:
      fMonTracksPerEvent = value;
      break;
    case GeantConfig::kMonTracks:
      fMonTracks = value;
    }
  }

  /** @brief Function returning the number of monitored features */
  inline int GetMonFeatures() const{
    // Get the number of monitored features
    return (fMonQueue + fMonMemory + fMonBasketsPerVol + fMonVectors + fMonConcurrency + fMonTracksPerEvent + fMonTracks);
  }

};
} // GEANT_IMPL_NAMESPACE
} // Geant
#endif // GEANT_RUN_CONFIG_H
