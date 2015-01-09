//===--- GeantScheduler.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantScheduler.h
 * @brief Implementation of dispatcher running in a single thread. Collects tracks
 * from all threads via an input queue and fills baskets corresponding to each
 * volume, which are then injected in the main work queue.
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_SCHEDULER
#define GEANT_SCHEDULER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifdef __STAT_DEBUG
#include "GeantTrackStat.h"
#endif

#if __cplusplus >= 201103L
#include <atomic>
#endif

#include "TMutex.h"

class concurrent_queue;
class GeantTrack;
class GeantBasket;
class GeantBasketMgr;

/**
 * @brief Class GeantScheduler
 * @details Dispatcher running in a single thread. Collects tracks
 * from all threads via an input queue and fills baskets corresponding to each
 * volume, which are then injected in the main work queue.
 * 
 */
class GeantScheduler : public TObject {
protected:
  Int_t fNvolumes;                   /** Number of active volumes in the geometry */
  Int_t fNpriority;                  /** Number of priority baskets held */
  GeantBasketMgr **fBasketMgr;       /** Array of basket managers */
  GeantBasketMgr *fGarbageCollector; /** Garbage collector manager */
#if __cplusplus >= 201103L
  std::atomic_int *fNtracks; /**[fNvolumes] Number of tracks per volume */
#endif
  Int_t fPriorityRange[2]; /** Prioritized events */
#ifdef __STAT_DEBUG
  GeantTrackStat fPStat; /** Statistics for the pending tracks */
  GeantTrackStat fQStat; /** Statistics for the queued tracks */
  GeantTrackStat fTStat; /** Statistics for the transported tracks */
#endif

private:

  /** @brief GeantScheduler constructor */
  GeantScheduler(const GeantScheduler &);

  /** @brief Function operator= */
  GeantScheduler &operator=(const GeantScheduler &);

public:

  /** GeantScheduler constructor */
  GeantScheduler();

  /** GeantScheduler destructor */
  virtual ~GeantScheduler();

  /**
   * @brief Function that add track
   * 
   * @param track Track that should be added
   */
  Int_t AddTrack(GeantTrack &track);

  /**
   * @brief Function for addition tracks
   * 
   * @param output Output GeantBasket object
   * @param ntot Total number of tracks
   * @param nnew Number of new tracks
   * @param nkilled Number of killed tracks
   */
  Int_t AddTracks(GeantBasket *output, Int_t &ntot, Int_t &nnew, Int_t &nkilled);

  /** @brief Function that adjust basket size */
  void AdjustBasketSize();

  /** @brief Function that create baskets */
  void CreateBaskets();

  /** @brief Function that clean baskets */
  void CleanBaskets();

  /** @brief Function that provides collection prioritized tracks*/
  Int_t CollectPrioritizedTracks();

  /**
   * @brief Function that return basket managers
   * return Array of basket managers
   * */
  GeantBasketMgr **GetBasketManagers() const { return fBasketMgr; }

  /**
   * @brief Function that return garbage collector
   * @return Garbage collector manager
   */
  GeantBasketMgr *GetGarbageCollector() const { return fGarbageCollector; }

#if __cplusplus >= 201103L

  /**
   * @brief Function that return N tracks
   * 
   * @param ib ?????
   * @return return number of tracks per volume
   */
  Int_t GetNtracks(Int_t ib) { return fNtracks[ib].load(); }
#endif

  /**
   * @brief Function that return N priority tracks
   * @return Number of priority baskets held
   */
  Int_t GetNpriority() const { return fNpriority; }

  /**
   * @brief Function that return N volumes
   * @return Number of active volumes in the geometry
   */
  Int_t GetNvolumes() const { return fNvolumes; }

  /**
   * @brief Function that return of range of prioritized events
   * 
   * @param min Minimum value of range
   * @param max Maximum value of range
   */
  void SetPriorityRange(Int_t min, Int_t max) {
    fPriorityRange[0] = min;
    fPriorityRange[1] = max;
  }
#ifdef __STAT_DEBUG

  /**
   * @brief [Function that return statistic of pending tracks
   * @return Statistics for the pending tracks
   */
  GeantTrackStat &GetPendingStat() { return fPStat; }

  /**
   * @brief Function that return statistic of queued tracks
   * @return Statistics for the queued tracks
   */
  GeantTrackStat &GetQueuedStat() { return fQStat; }

  /**
   * @brief Function that return statistic of transported tracks
   * @return Statistics for the transported tracks
   */
  GeantTrackStat &GetTransportStat() { return fTStat; }
#endif

  /** @brief Function that flush list of priority tracks */
  Int_t FlushPriorityBaskets();

  /** @brief Garbage collection function */
  Int_t GarbageCollect();

  /** @brief Function that print size */
  void PrintSize() const;

  /** @brief Function that returns size */
  size_t Sizeof() const;

  ClassDef(GeantScheduler, 1) // Main basket scheduler
};
#endif
