//===--- GeantBasket.h - Geant-V --------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantBasket.h
 * @brief Implementation of baskets of tracks for Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_BASKET
#define GEANT_BASKET

#ifdef USE_ROOT
#include "TGeoExtension.h"
#else
#define MIC_BIT(n) (1ULL<<(n))
#endif

#include "Geant/Typedefs.h"
#include "GeantTrackVec.h"
#include "priority_queue.h"
#include "mpmc_bounded_queue.h"
#include "GeantPropagator.h"
#include "GeantFwd.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantBasketMgr;

/**
 * @brief Class GeantBasket descripting basic operations with baskets
 * @details Basket of tracks in the same volume which are transported by a single thread
 */
class GeantBasket {
protected:
  GeantBasketMgr *fManager; /** Manager for the basket */
public:
  std::atomic_int fNcopying;    /** Number tracks copying concurrently */
  std::atomic_int fNbooked;     /** Number of slots booked for copying */
  std::atomic_int fNcopied;     /** Number of tracks copied */
  std::atomic_int fNused;       /** Number of threads using the basket */
  size_t fIbook0;               /** Start slot number */
  std::atomic_flag fDispatched; /** Atomic flag marking the basket as dispatched */
  int fThreshold;               /** Current transport threshold */
protected:
  // GeantHit_v        fHits;  /** Vector of produced hits */
  GeantTrack_v fTracksIn;  /** Vector of input tracks */
//  GeantTrack_v fTracksOut; /** Vector of output tracks */

private:
  /** @todo Still not implemented */
  GeantBasket(const GeantBasket &);

  /** @todo Still not implemented operator = */
  GeantBasket &operator=(const GeantBasket &);

  bool fIsMixed = false;

public:
  /** @brief Default GeantBasket constructor */
  GeantBasket();

  /**
   * @brief GeantBasket standard constructor
   *
   * @param size Initial size of input/output track arrays
   * @param mgr  Basket manager handling this basket
   */
  GeantBasket(GeantPropagator* prop, int size, GeantBasketMgr *mgr);

  /**
   * @brief GeantBasket standard constructor
   *
   * @param size Initial size of input/output track arrays
   * @param depth Maximum geometry depth
   */
  GeantBasket(int size, int depth);

  /** @brief GeantBasket destructor */
  virtual ~GeantBasket();

  /**
   * @brief Add a scalar track to basket input.
   * @details Add concurrently track from generator or physics process
   *
   * @param track Reference to track object
   * @return Track index;
   */
  int AddTrack(GeantTrack &track);

  /**
   * @bref Add a track from vector container to basket input.
   * @details Add concurrently a track to basket input.
   *
   * @param tracks Array of tracks to copy from.
   * @param itr Track id.
   * @return Track index;
   */
  int AddTrack(GeantTrack_v &tracks, int itr);

  /**
   * @brief Function to add multiple tracks to basket
   * @details Add multiple tracks from a track_v array
   *
   * @param tracks Tracks from a track_v array
   * @param istart Start track id
   * @param iend End track id
   */
  void AddTracks(GeantTrack_v &tracks, int istart, int iend);

  /** @brief Book a slot and return number of slots booked */
  int BookSlot() { return (fNbooked.fetch_add(1, std::memory_order_seq_cst) + 1); }

  /** @brief Get number of booked slots */
  int GetNbooked() { return (fNbooked.load()); }

  /** @brief Virtual function for clearing the basket */
  virtual void Clear(const char *option = "");

  /**
   * @brief Check if a basket contains tracks in a given event range
   *
   * @param evstart Start event id.
   * @param nevents Number of events (default 1)
   */
  bool Contains(int evstart, int nevents = 1) const;

  /**
   * @brief Function returning the number of input tracks
   * @return Number of input tracks
   */
  int GetNinput() const { return fTracksIn.GetNtracks(); }

  /**
   * @brief Function returning the number of output tracks
   * @return Number of output tracks
   */
//  int GetNoutput() const { return fTracksOut.GetNtracks(); }

  /**
   * @brief Function returning a reference to the vector of input tracks
   * @return Reference to input vector of tracks
   */
  GeantTrack_v &GetInputTracks() { return fTracksIn; }

  /**
   * @brief Function returning a reference to the vector of output tracks
   * @return Reference to output vector of tracks
   */
//  GeantTrack_v &GetOutputTracks() { return fTracksOut; }

  /**
   * @brief Function returning the manager of the basket
   * @return Pointer to manager of the basket
   */
  GeantBasketMgr *GetBasketMgr() const { return fManager; }

  /**
   * @brief Function setting the manager of the basket
   */
  void SetBasketMgr(GeantBasketMgr *mgr) { fManager = mgr; }

  /**
   * @brief Function for defining basket transportability threshold
   * @return  Value of transportability threshold
   */
  int GetThreshold() const { return fThreshold; }

  /**
   * @brief Function returning the volume for this basket
   * @return Pointer to associated logical volume
   */
  Volume_t *GetVolume() const;

  /**
   * @brief Function returning the mixed tracks property.
   * @return Boolean value if the tracks are mixed from several volumes
   */
  bool IsMixed() const { return fIsMixed; }
  /**
   * @brief Check if tracks are being copied
   * @return Track copy flag
   */
  inline bool IsCopying() const { return (fNcopying.load(std::memory_order_seq_cst)); }
  //  inline bool IsCopying() const { return fNcopying; }
  inline int GetNcopying() const { return fNcopying.load(); }
  //  inline int GetNcopying() const { return fNcopying; }

  /**
   * @brief Mark start of copy operation
   * @return Number of concurrent track copy operations
   */
  inline int StartCopying() { return (fNcopying.fetch_add(1, std::memory_order_seq_cst) + 1); }
  //  inline int StartCopying() { return ( ++fNcopying ); }

  /**
   * @brief Mark stop of copy operation
   * @return Number of concurrent track copy operations remaining
   */
  inline int StopCopying() { return (fNcopying.fetch_sub(1, std::memory_order_seq_cst) - 1); }
  //  inline int StopCopying() { return ( --fNcopying ); }

  /**
   * @brief Print the basket content
   */
  virtual void Print(const char *option = "") const;

  /**
   * @brief Print the parameters for a given track
   *
   * @param itr Track id.
   * @param input Refer to input or output track (default input)
   */
  void PrintTrack(int itr, bool input = true) const;

  /** @brief Recycle this basket */
  void Recycle(GeantTaskData *td);

  /**
   * @brief  Function that providing the size of this basket in bytes
   * @return Sum of sizes of all tracks (input and output) + data members
   */
  size_t Sizeof() const {
    return fTracksIn.Sizeof() + sizeof(this) + sizeof(GeantBasketMgr *) +
           sizeof(std::atomic_int);
  }

  /**
   * @brief Flag the basket to contain tracks that are mixed from different volumes
   *
   * @param flag Boolean flag.
   */
  void SetMixed(bool flag) { fIsMixed = flag; }

  /**
   * @brief Function to change the transportability threshold for the basket
   *
   * @param threshold New threshold value
   */
  void SetThreshold(int threshold);

  /** @brief Try to get the dispatch lock for the basket */
  bool TryDispatch() { return (!fDispatched.test_and_set(std::memory_order_acquire)); }

};

class GeantScheduler;

/**
 * @brief Class managing all baskets for a given logical volume
 * @details Basket manager for a given volume. Holds a list of free baskets stored in a
 * concurrent queue, a current basket to be filled and a priority basket used
 * in priority mode.
 */
#ifdef USE_ROOT
class GeantBasketMgr : public TGeoExtension {
#else
class GeantBasketMgr {
#endif

public:
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;

protected:
  GeantScheduler *fScheduler; /** Scheduler for this basket */
  Volume_t *fVolume;          /** Volume for which applies */
  int fNumber;                /** Number matching the volume index */
  int fBcap;                  /** Maximum capacity of baskets held */
  int fQcap;                  /** Queue capacity */
  bool fActive;               /** Activity flag for generating baskets */
  bool fCollector;            /** Mark this manager as event collector */
  std::atomic_int fThreshold; /** Adjustable transportability threshold */
  std::atomic_int fNbaskets;  /** Number of baskets for this volume */
  std::atomic_int fNused;     /** Number of baskets in use */
  std::atomic<size_t> fIbook; /** Booking index */
  typedef std::atomic<GeantBasket *> atomic_basket;
  atomic_basket fCBasket;                        /** Current basket being filled */
  Geant::priority_queue<GeantBasket *> *fFeeder; /** Feeder queue to which baskets get injected */
  std::vector<GeantBasket *> fDispatchList;      /** list of baskets to be dispatched */

private:
  /** @todo Still not implemented */
  GeantBasketMgr(const GeantBasketMgr &);

  /** @todo Still not implemented operator = */
  GeantBasketMgr &operator=(const GeantBasketMgr &);

  /**
   * @brief The caller thread books the basket held by the atomic for track addition
   *
   * @param current Current atomic basket to be booked
   * @return Number of booked slots
   * @return Basket actually booked
   */
  GeantBasket *BookBasket(GeantTaskData *td);

  /**
   * @brief Attempt to replace the atomic content of the basket with new one, strong CAS
   *
   * @param expected expected value for book index
   * @param td Task data
   * @return True if replacement made by the call
   */
  bool ReplaceBasketStrong(size_t expected, GeantTaskData *td);

  /**
   * @brief Attempt to replace the atomic content of the basket with new one, weak CAS
   *
   * @param expected expected value for book index
   * @param td Task data
   * @return True if replacement made by the call
   */
  bool ReplaceBasketWeak(size_t expected, GeantTaskData *td);

public:
  /**
   * @brief Get next free basket
   * @details Get pointer to next free basket
   */
  GeantBasket *GetNextBasket(GeantTaskData *td);

  /**
   * @brief Create a number of empty baskets
   */
  void CreateEmptyBaskets(int nbaskets, GeantTaskData *td);

public:
  /** @brief GeantBasketMgr dummy constructor */
  GeantBasketMgr()
      : fScheduler(0), fVolume(0), fNumber(0), fBcap(0), fQcap(0), fActive(false), fCollector(false), fThreshold(0),
        fNbaskets(0), fNused(0), fIbook(0), fCBasket(0), fFeeder(0), fDispatchList() {}

  /** @brief GeantBasketMgr normal constructor
   *
   * @param sch Scheduler dealing with this basket manager
   * @param vol Volume associated with this
   * @param number Number for the basket manager
  */
  GeantBasketMgr(GeantPropagator* prop, GeantScheduler *sch, Volume_t *vol, int number, bool collector = false);

  /** @brief Destructor of GeantBasketMgr */
  virtual ~GeantBasketMgr();

  /**
   * @brief Grab function
   * @details Interface of TGeoExtension for getting a reference to this from Volume
   * @return Pointer to the base class
   */
 #ifdef USE_ROOT
  virtual TGeoExtension *Grab() { return this; }
 #endif 

  /**
   * @brief Release function
   * @details Interface of TGeoExtension to signal releasing ownership of this from TGeoVolume
   */
  virtual void Release() const {}

  /**
   * @brief Activate this basket manager for generating baskets
   * @details Activation happens on threshold on percentage of tracks transported
   * in the associated volume and threshold on total number of tracks (learning phase)
   */
  void Activate(GeantPropagator* prop);

  /**
   * @brief Function adding a track to basket up to the basket threshold
   *
   * @param track  Track that should be added to basket
   * @param priority Set priority (by default false)
   */
  int AddTrack(GeantTrack &track, bool priority, GeantTaskData *td);

  /**
   * @brief Function adding an indexed track to basket up to the basket threshold
   *
   * @param trackv Array of tracks containing the track to be added
   * @param itr Track id
   * @param priority Set priority (by default false)
   */
  int AddTrack(GeantTrack_v &trackv, int itr, bool priority, GeantTaskData *td);

  /**
   * @brief Function adding by a unique thread a track to the basket manager up to the basket threshold
   *
   * @param trackv Array of tracks containing the track to be added
   * @param itr Track id
   * @param priority Set priority (by default false)
   */
  int AddTrackSingleThread(GeantTrack_v &trackv, int itr, bool priority, GeantTaskData *td);

  /**
   * @brief Function cleaning a number of free baskets
   *
   * @param ntoclean Number of baskets to be cleaned
   */
  void CleanBaskets(int ntoclean, GeantTaskData *td);

  /** @brief Function doing full collection to work queue of non-empty baskets*/
  int GarbageCollect(GeantTaskData *td);

  /**
   * @brief Function that set capacity of baskets
   *
   * @param capacity Capacity of baskets to be set
   */
  void SetBcap(int capacity) { fBcap = capacity; }

  /**
   * @brief Function that returns the capacity of baskets
   * @return Maximum capacity of baskets held
   */
  int GetBcap() const { return fBcap; }

  /**
   * @brief Snapshot of the number of baskets
   * @return Number of baskets
   */
  int GetNbaskets() const { return fNbaskets.load(); }

  /**
   * @brief Snapshot of the number of baskets in use
   * @return number of baskets in use
   */
  int GetNused() const { return fNused.load(); }

  /**
   * @brief Snapshot of the current basket threshold
   * @return Threshold value
   */
  int GetThreshold() const { return fThreshold.load(); }

  /**
   * @brief Check if the basket manager has tracks filled and pending threshold
   * @return True if the current basket holds tracks
   */
  bool HasTracks() const { return GetCBasket()->GetNinput(); }

  /**
   * @brief Getter for activity (manager generating baskets)
   * @return Value for activity flag
   */
  bool IsActive() const { return fActive; }

  /**
   * @brief Getter for garbage collector status
   * @return Value for collector flag
   */
  bool IsCollector() const { return fCollector; }

  /** @brief Setter for garbage collector status */
  void SetCollector() { fCollector = true; }

  /**
   * @brief Function to set transportability threshold
   *
   * @param thr Threshold that should be set
   */
  void SetThreshold(int thr) { fThreshold.store(thr); }

  /**
   * @brief Function that will load the current basket
   * @return Load current basket
   */
  GeantBasket *GetCBasket() const { return fCBasket.load(); }

  /**
   * @brief Get the current basket in non-thread safe way to send for transport.
   * @details Method called only by the thread owning this manager that should
   * have fCollector=true
   * @return Current basket
   */
  GeantBasket *GetBasketForTransport(GeantTaskData *td) {
    GeantBasket *basket = GetCBasket();
    SetCBasket(GetNextBasket(td));
    return basket;
  }

  /**
   * @brief Function to set current basket
   *
   * @param basket Basket to be set as current
   */
  void SetCBasket(GeantBasket *basket) { fCBasket.store(basket); }

  /**
   * @brief Function that returns the scheduler
   * @return Scheduler for basket
   */
  GeantScheduler *GetScheduler() const { return fScheduler; }

  /** @brief Function that returns the name of volume */
  #ifdef USE_ROOT
  const char *GetName() const { return (fVolume) ? fVolume->GetName() : ClassName(); }
  #else
  const char *GetName() const { return (fVolume) ? fVolume->GetName() : "No Volume"; }
  #endif
  /**
   * @brief Function that returns the number assigned to basket
   * @return Number assigned to basket
   */
  int GetNumber() const { return fNumber; }

  /**
   * @brief Function that returns the associated volume pointer
   * @return Volume for which applies basket
   */
  Volume_t *GetVolume() const { return fVolume; }

  /** @brief Print the current basket */
  virtual void Print(const char *option = "") const;

  /**
   * @brief Push the basket to the queue recalculating threshold
   * @param basket Basket to be pushed
   * @param priority Priority on or off
   */
  void Push(GeantBasket *basket, bool priority, GeantTaskData *td);

  /**
   * @brief Recycles a given basket
   *
   * @param b Pointer to current GeantBasket for recycling
   */
  void RecycleBasket(GeantBasket *b, GeantTaskData *td);

  /**
   * @brief Function setting the feeder work queue
   *
   * @param queue priority_queue for GeantBasket
   */
  void SetFeederQueue(Geant::priority_queue<GeantBasket *> *queue) { fFeeder = queue; }

  /**
   * @brief Function returning the size of the basket being filled
   * @return Size of basket
   */
  size_t Sizeof() const {
    GeantBasket *c = GetCBasket();
    return (c) ? (sizeof(GeantBasketMgr) + (GetNbaskets()) * c->Sizeof()) : sizeof(GeantBasketMgr);
  }

  /** @brief Print size */
  void PrintSize() const;

  /**
   * @brief Function GetFeederQueue()
   * @return Feeder queue to which baskets get injected
   */
  Geant::priority_queue<GeantBasket *> *GetFeederQueue() const { return fFeeder; }

 #ifdef USE_ROOT
  ClassDef(GeantBasketMgr, 0) // A path in geometry represented by the array of indices
 #endif 
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
