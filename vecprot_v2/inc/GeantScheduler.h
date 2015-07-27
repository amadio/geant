//===--- GeantScheduler.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantScheduler.h
 * @brief Implementation of the GeantV scheduler running in a single thread. Collects tracks
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

#include <atomic>
#include "TMutex.h"

class concurrent_queue;
class GeantBasket;
class GeantBasketMgr;

#include "GeantFwd.h"

/**
 * @brief Class GeantScheduler
 * @details Dispatcher running in a single thread. Collects tracks
 * from all threads via an input queue and fills baskets corresponding to each
 * volume, which are then injected in the main work queue.
 * 
 */
class GeantScheduler : public TObject {
public:
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;
protected:
  int fNvolumes;                   /** Number of active volumes in the geometry */
  int fNpriority;                  /** Number of priority baskets held */
  GeantBasketMgr **fBasketMgr;       /** Array of basket managers */
  GeantBasketMgr *fGarbageCollector; /** Garbage collector manager */
  int           *fNstvol;  /**[fNvolumes] Number of steps per volume */
  int           *fIstvol;  /**[fNvolumes] Sorted index of number of steps per volume */
  int           *fNvect;   /**[256] Number of tracks basketized in vectors of given size */
  std::atomic_int  fNsteps;  /** Total number of tracks steps */
  std::atomic_int  fCrtMgr;  /** Current basket manager being garbage collected */
  std::atomic_bool fCollecting;      /** Flag marking colecting tracks for priority events */
  std::atomic_flag fLearning;        /** Flag marking the learning phase */
  std::atomic_flag fGBCLock;         /** Flag marking that garbage collector is busy */
  int fPriorityRange[2]; /** Prioritized events */

private:

  /** 
   * @brief GeantScheduler copy constructor 
   * @details Not allowed
   */
  GeantScheduler(const GeantScheduler &);

  /** 
   * @brief Function operator= 
   * @details Not allowed
   */
  GeantScheduler &operator=(const GeantScheduler &);

public:

  /** GeantScheduler default constructor */
  GeantScheduler();

  /** GeantScheduler destructor */
  virtual ~GeantScheduler();

  /** @brief Activate basket managers based on distribution of steps */
  void ActivateBasketManagers();

  /**
   * @brief Schedule a new track
   * 
   * @param track Track to be scheduled
   */
  int AddTrack(GeantTrack &track, GeantTaskData *td);

  /**
   * @brief Re-schedule all tracks from an output basket
   * @details Transported baskets contain tracks exiting the current volume,
   * tracks killed and new tracks generated by physics along the step. These tracks
   * have to be re-scheduled for transport or bookkeeped for removal. 
   *
   * @param output Transported basket
   * @param ntot Total number of tracks
   * @param nnew Number of new tracks
   * @param nkilled Number of killed tracks
   * @param td Thread data
   */
  int AddTracks(GeantBasket *output, int &ntot, int &nnew, int &nkilled, GeantTaskData *td);

  /** @brief Function to adjust the basket size automatically */
  void AdjustBasketSize();

  /** @brief Function to create initially baskets */
  void CreateBaskets();

  /**
   * @brief Getter for the array of basket managers
   * return Array of basket managers
   * */
  GeantBasketMgr **GetBasketManagers() const { return fBasketMgr; }

  /**
   * @brief Getter for the garbage collector
   * @return Garbage collector manager
   */
  GeantBasketMgr *GetGarbageCollector() const { return fGarbageCollector; }

  /**
   * @brief Getter for total number of steps
   * 
   * @return Number of steps
   */
  int GetNsteps() const { return fNsteps.load(); }

  /**
   * @brief Getter for collecting flag
   * @return Value of fCollecting flag
   */
   Bool_t IsCollecting() const { return fCollecting.load(); }
   
  /**
   * @brief Setter for collecting flag
   */
   void SetCollecting(Bool_t flag) { fCollecting.store(flag); }

  /**
   * @brief Getter for learning flag
   * @return Value of fLearning flag
   */
   Bool_t IsLearning() { 
     bool learning = fLearning.test_and_set(std::memory_order_acquire);
     if (!learning) fLearning.clear(std::memory_order_release);
     return learning; }
   
  /**
   * @brief Setter for the learning flag
   */
   void SetLearning(Bool_t flag) { 
     if (flag) fLearning.test_and_set(std::memory_order_acquire);
     else fLearning.clear(std::memory_order_release); }

  /**
   * @brief Getter for array fNvect
   * @return Pointer to fNvect array
   */
  int *GetNvect() { return fNvect; } 
   
  /**
   * @brief Function to return N priority baskets per volume
   * @return Number of priority baskets held
   */
  int GetNpriority() const { return fNpriority; }

  /**
   * @brief Function to return N volumes
   * @return Number of active volumes in the geometry
   */
  int GetNvolumes() const { return fNvolumes; }

  /**
   * @brief Function to return the range of prioritized events
   * 
   * @param min Minimum value of range
   * @param max Maximum value of range
   */
  void SetPriorityRange(int min, int max) {
    fPriorityRange[0] = min;
    fPriorityRange[1] = max;
  }

  /** @brief Garbage collection function */
  int GarbageCollect(GeantTaskData *td, Bool_t force=false);

  /** @brief Function to print size */
  void PrintSize() const;

  /** @brief Function to returns size */
  size_t Sizeof() const;

  ClassDef(GeantScheduler, 1) // Main basket scheduler
};
#endif
