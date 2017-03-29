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

  int fNtotal;            /** Total number of events to be transported */
  int fNbuff;             /** Number of buffered events */
  int fNprocesses;        /** Number of active physics processes */
  int fNstart;            /** Cumulated initial number of tracks */
  int fMaxTracks;         /** Maximum number of tracks per event */
  int fMaxThreads;        /** Maximum number of threads */
  int fNminThreshold;     /** Threshold for starting transporting a basket */
  int fNvecThreshold;     /** Threshold for executing in vector mode */
  int fDebugEvt;          /** Event to debug */
  int fDebugTrk;          /** Track to debug */
  int fDebugStp;          /** Step to start debugging */
  int fDebugRep;          /** Number of steps to debug */
  int fMaxSteps;          /** Maximum number of steps per track */
  int fNperBasket;        /** Number of tracks per basket */
  int fMaxPerBasket;      /** Maximum number of tracks per basket */
  int fMaxPerEvent;       /** Maximum number of tracks per event */
  int fMaxDepth;          /** Maximum geometry depth */
  int fLearnSteps;        /** Number of steps needed for the learning phase */
  int fLastEvent;         /** Last transported event */
  float fPriorityThr;     /** Threshold for prioritizing events */
  int fNstepsKillThr;     /** Threshold in number of steps to kill a track */
  int fNminReuse;         /** Minimum number of transported tracks to be reused without re-basketizing */
  int fNstackLanes;       /** Number of stacked lanes in the stack-like buffers */
  int fNmaxBuffSpill;     /** Maximum number of tracks spilled from the stack-like buffer per stepping iteration */

  double fMaxRes;         /** Maximum resident memory allowed [MBytes] */
  double fMaxVirt;        /** Maximum virtual memory allowed [MBytes] */
  double fNaverage;       /** Average number of tracks per event */
  double fVertex[3];      /** Vertex position */
  double fEmin;           /** Min energy threshold */
  double fEmax;           /** Max energy threshold */
  double fBmag;           /** Magnetic field */
  double fEpsilonRK;      /** Relative error in RK integration */

  bool fUsePhysics;       /** Enable/disable physics */
  bool fUseRungeKutta;    /** Enable/disable Runge-Kutta integration in field */
  bool fUseDebug;         /** Use debug mode */
  bool fUseGraphics;      /** Graphics mode */
  bool fUseStdScoring;    /** Use standard scoring */
  bool fUseV3;            /** Use version 3 of the scheduler */

  int fMonQueue;          /** Monitor the work queue */
  int fMonMemory;         /** Monitor the memory */
  int fMonBasketsPerVol;  /** Monitor baskets per volume */
  int fMonVectors;        /** Monitor vector scheduling */
  int fMonConcurrency;    /** Monitor concurrency */
  int fMonTracksPerEvent; /** Monitor tracks status per event */
  int fMonTracks;         /** Monitor number of tracks */

  bool fFillTree;         /** Enable I/O */
  bool fUseMonitoring;    /** Monitoring different features */
  bool fUseAppMonitoring; /** Monitoring the application */
  int  fTreeSizeWriteThreshold; /** Maximum size of the tree (before automatic writing) **/
  bool fConcurrentWrite;  /** Switch between single and mutlithreaded writing */
#ifndef VECCORE_CUDA
  std::string fGeomFileName; /** Geometry file name */
#endif

public:
  VECCORE_ATT_DEVICE
  GeantConfig(): fNtotal(1000), fNbuff(100), fNprocesses(3), fNstart(0), fMaxTracks(0), fMaxThreads(100), fNminThreshold(16), fNvecThreshold(8),
    fDebugEvt(-1), fDebugTrk(-1), fDebugStp(-1), fDebugRep(-1), fMaxSteps(10000), fNperBasket(16), fMaxPerBasket(256),
    fMaxPerEvent(0), fMaxDepth(0), fLearnSteps(0), fLastEvent(0), fPriorityThr(0), fNstepsKillThr(50000),
    fNminReuse(10000), fNstackLanes(10), fNmaxBuffSpill(128), fMaxRes(0), fMaxVirt(0), fNaverage(0), fVertex(),
    fEmin(1.E-4),// 100 KeV
    fEmax(10),// 10 Gev
    fBmag(0.),// kiloGauss
    fEpsilonRK(0.0003), fUsePhysics(true), fUseRungeKutta(false), fUseDebug(false), fUseGraphics(false), fUseStdScoring(false), fUseV3(false),
    fMonQueue(0), fMonMemory(0), fMonBasketsPerVol(0), fMonVectors(0), fMonConcurrency(0),
    fMonTracksPerEvent(0), fMonTracks(0), fFillTree(false), fUseMonitoring(false), fUseAppMonitoring(false),
    fTreeSizeWriteThreshold(100000), fConcurrentWrite(true)
#ifndef VECCORE_CUDA
      , fGeomFileName()
#endif
                 {};

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
