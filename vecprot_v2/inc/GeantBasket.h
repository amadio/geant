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

#include "TObject.h"
#include "TGeoExtension.h"
#include "TGeoVolume.h"
#include "GeantTrack.h"
#include "priority_queue.h"
#include "mpmc_bounded_queue.h"

class TGeoVolume;
class GeantBasketMgr;

/**
 * @brief Class GeantBasket descripting basic operations with baskets
 * @details Basket of tracks in the same volume which are transported by a single thread
 */
class GeantBasket : public TObject {
public:

  /**
   * @enum EbasketFlags
   * @details Basked mixing tracks from different volumes
   */
  enum EBasketFlags {
    kMixed = BIT(14) 
  };

protected:
  GeantBasketMgr *fManager; /** Manager for the basket */
  GeantTrack_v fTracksIn;   /** Vector of input tracks */
  GeantTrack_v fTracksOut;  /** Vector of output tracks */
//GeantHit_v        fHits;  /** Vector of produced hits */
#if __cplusplus >= 201103L
  std::atomic_int fAddingOp; /** Number of track adding ops */
#endif
  Int_t fThreshold; /** Transport threshold */

private:

  /** @todo Still not implemented */
  GeantBasket(const GeantBasket &);

  /** @todo Still not implemented operator = */
  GeantBasket &operator=(const GeantBasket &);
public:
  
  /** @brief Simple GeantBasket constructor */
  GeantBasket();
  
  /**
   * @brief GeantBasket parameterized constructor 
   * 
   * @param size Size of created GeantBasket
   * @param mgr  GeantBasket manager for basket
   */
  GeantBasket(Int_t size, GeantBasketMgr *mgr);
  
  /** @brief GeantBasket destructor */
  virtual ~GeantBasket();

  /**
   * @brief Function that add track to basket
   * @details Add track from generator or physics process
   * 
   * @param track Track from generator or physics process
   */
  void AddTrack(GeantTrack &track);
  
  /**
   * @brief Function that add track to basket
   * @details Add track from a track_v array (copied)
   * 
   * @param tracks Track from a track_v array 
   * @param itr Track ID
   */
  void AddTrack(GeantTrack_v &tracks, Int_t itr);
  
  // Add multiple tracks from a track_v array
  /**
   * @brief Function that add multiple tracks to basket
   * @details Add multiple tracks from a track_v array
   * 
   * @param tracks Tracks from a track_v array 
   * @param istart Start track ID 
   * @param iend End track ID
   */
  void AddTracks(GeantTrack_v &tracks, Int_t istart, Int_t iend);
  
  /** @brief Virtual function for cleaning */
  virtual void Clear(Option_t *option = "");
  
  /**
   * @brief Function of basket containing
   * 
   * @param evstart  Start event ID ?????
   * @param nevents Number of events (by default 1)
   */
  Bool_t Contains(Int_t evstart, Int_t nevents = 1) const;
  
  /**
   * @brief Function that return number of N tracks in vector of input tracks
   * @return Number of N input tracks
   */
  Int_t GetNinput() const { return fTracksIn.GetNtracks(); }
  
  /**
   * @brief Function that return number of N tracks in vector of output tracks
   * @return Number of N output tracks
   */
  Int_t GetNoutput() const { return fTracksOut.GetNtracks(); }
  
  /**
   * @brief Function that return vector of input tracks
   * @return Vector of N input tracks
   */
  GeantTrack_v &GetInputTracks() { return fTracksIn; }
  
  /**
   * @brief Function for definition vector of output tracks
   * @return Vector of N output tracks
   */
  GeantTrack_v &GetOutputTracks() { return fTracksOut; }
  
  /**
   * @brief Function for definition basket manager 
   * @return Manager for the basket
   */
  GeantBasketMgr *GetBasketMgr() const { return fManager; }
  
  /**
   * @brief Function for definition threshold
   * @return  Value of transport threshold
   */
  Int_t GetThreshold() const { return fThreshold; }
  
  /**
   * @brief Function that return volume
   * @return TGeoVolume object
   */
  TGeoVolume *GetVolume() const;
  
  /**
   * @brief Function that return flag if tracks are mixed or not
   * @return Boolean value if ((fBits & kMixed) != 0);
   */
  Bool_t IsMixed() const { return TObject::TestBit(kMixed); }
  
  /**
   * @brief Load number of track adding ops
   * @return fAddingOp -> number of track adding ops
   */
  inline Bool_t IsAddingOp() const { return (fAddingOp.load()); }
  
  /**
   * @brief Lock number of track adding ops
   * @return fAddingOp -> Number of track adding ops
   */
  inline Int_t LockAddingOp() { return ++fAddingOp; }
  
  /**
   * @brief Unlock number of track adding ops
   * @return fAddingOp -> Number of track adding ops
   */
  inline Int_t UnLockAddingOp() { return --fAddingOp; }
  
  /**
   * @brief Simple print function
   */
  virtual void Print(Option_t *option = "") const;
  
  /**
   * @brief Function of printing track
   * 
   * @param itr Track ID
   * @param input Flag that checks if input exist
   */
  void PrintTrack(Int_t itr, Bool_t input = kTRUE) const;
  
  /** @brief Recycle function */
  void Recycle();
  
  /**
   * @brief  Function that provides size of basket
   * @return Sum of tracks inside of basket & outside tracks & size of GeantBasketMgr & size of atomic<int>
   */
  size_t Sizeof() const {
    return fTracksIn.Sizeof() + fTracksOut.Sizeof() + sizeof(TObject) + sizeof(GeantBasketMgr *) +
           sizeof(std::atomic_int);
  }

  /**
   * @brief Function that set flag if tracks are mixed
   * 
   * @param flag Flag that provides kMixed 
   */ 
  void SetMixed(Bool_t flag) { TObject::SetBit(kMixed, flag); }
  
  /**
   * @brief Function that set threshold for basket
   * 
   * @param threshold Threshold
   */
  void SetThreshold(Int_t threshold);

  ClassDef(GeantBasket, 1) // A basket containing tracks in the same geomety volume
};

class GeantScheduler;

/**
 * @brief Class of basket manager
 * @details Basket manager for a given volume. Holds a list of free baskets stored in a
 * concurrent queue
 */
class GeantBasketMgr : public TGeoExtension {
protected:
  GeantScheduler *fScheduler; /** Scheduler for this basket */
  TGeoVolume *fVolume;        /** Volume for which applies */
  Int_t fNumber;              /** Number assigned */
  Int_t fBcap;                /** Max capacity of baskets held */
  Int_t fQcap;                /** Queue capacity */
#if __cplusplus >= 201103L
  std::atomic_int fThreshold; /** Adjustable transportability threshold */
  std::atomic_int fNbaskets;  /** Number of baskets for this volume */
  std::atomic_int fNused;     /** Number of baskets in use */
  typedef std::atomic<GeantBasket *> atomic_basket;
  atomic_basket fCBasket;  /** Current basket being filled */
  atomic_basket fPBasket;  /** Current priority basket being filled */
  std::atomic_flag fLock;  /** Atomic lock for stealing current basket */
  std::atomic_flag fQLock; /** Atomic lock for increasing queue size */
#endif
  mpmc_bounded_queue<GeantBasket *> *fBaskets;   /** Queue of available baskets */
  Geant::priority_queue<GeantBasket *> *fFeeder; /** Feeder queue to which baskets get injected */
  TMutex fMutex;                                 /** Mutex for this basket manager */
private:

  /** @todo Still not implemented */
  GeantBasketMgr(const GeantBasketMgr &);

  /** @todo Still not implemented operator = */          
  GeantBasketMgr &operator=(const GeantBasketMgr &);
#if __cplusplus >= 201103L

  /**
   * @brief Steal and Replace function
   *  
   * @param current Current atomic basket
   */
  GeantBasket *StealAndReplace(atomic_basket &current);

  /**
   * @brief Steal and Pin function
   * 
   * @param current Current atomic basket
   */
  GeantBasket *StealAndPin(atomic_basket &current);

  /**
   * @brief Steal matching function
   * 
   * @param global Global atomic basket
   * @param content Content of GeantBasket 
   */
  Bool_t StealMatching(atomic_basket &global, GeantBasket *content);
#endif

  /**
   * @brief Function increasing queue size 
   * 
   * @param newsize New size of queue  
   */
  Bool_t IncreaseQueueSize(Int_t newsize);

public:

  /**
   * @brief Get next basket 
   * @details Get pointer to next Basket
   */
  GeantBasket *GetNextBasket();

public:

  /** @brief GeantBasketMgr constructor */
  GeantBasketMgr()
      : fScheduler(0), fVolume(0), fNumber(0), fBcap(0), fQcap(0), fThreshold(0), fNbaskets(0),
        fNused(0), fCBasket(0), fPBasket(0), fLock(), fQLock(), fBaskets(0), fFeeder(0), fMutex() {}
  GeantBasketMgr(GeantScheduler *sch, TGeoVolume *vol, Int_t number);
  
  /** @brief Destructor of GeantBasketMgr */
  virtual ~GeantBasketMgr();
  
  /**
   * @brief Grab function
   * @details Function to interface for getting a reference (grab)
   * @return Pointer to the extension
   */
  virtual TGeoExtension *Grab() { return this; }

  /** @brief Method called always when the pointer to the extension is not needed */
  virtual void Release() const {}

  /**
   * @brief Function that add track to basket
   * 
   * @param track  Track that should be added to basket
   * @param priority Set priority (by default kFALSE)
   */
  Int_t AddTrack(GeantTrack &track, Bool_t priority = kFALSE);

  /**
   * @brief Function that add track to basket
   * 
   * @param trackv Track from track_v array that should be added to basket
   * @param itr Track ID
   * @param priority Set priority (by default kFALSE)
   */
  Int_t AddTrack(GeantTrack_v &trackv, Int_t itr, Bool_t priority = kFALSE);

  /**
   * @brief Collection of prioterized tracks
   * 
   * @param evmin Minimum quantity of events
   * @param evmax Maximum quantity of events
   */
  Int_t CollectPrioritizedTracks(Int_t evmin, Int_t evmax);

  /**
   * @brief Function that clean baskets
   * 
   * @param ntoclean Number of baskets to be cleaned
   */
  void CleanBaskets(Int_t ntoclean);

  /** @brief Function that flush only priority baskets */
  Int_t FlushPriorityBasket();

  /** @brief Function of collection of garbage */
  Int_t GarbageCollect();

  /**
   * @brief Function that set capacity of baskets
   * 
   * @param capacity Capacity of baskets tjat should be set
   */
  void SetBcap(Int_t capacity) { fBcap = capacity; }

  /**
   * @brief Function that return capacity of baskets
   * @return Maximum capacity of baskets held
   */
  Int_t GetBcap() const { return fBcap; }

#if __cplusplus >= 201103L

   /**
    * @brief Get number of baskets
    * @return Number of baskets
    */
  Int_t GetNbaskets() const { return fNbaskets.load(); }

  /**
   * @brief Get number of baskets in use
   * @return number of baskets in use
   */
  Int_t GetNused() const { return fNused.load(); }

  /**
   * @brief Function that get threshold 
   * @return Threshold that should be load
   */
  Int_t GetThreshold() const { return fThreshold.load(); }

  /**
   * @brief Function that set threshold
   * 
   * @param thr Threshold that should be set
   */
  void SetThreshold(Int_t thr) { fThreshold.store(thr); }

  /**
   * @brief Function that will load current basket
   * @return Load current basket
   */
  GeantBasket *GetCBasket() const { return fCBasket.load(); }

  /**
   * @brief Function that will load priority basket
   * @return Load priority basket 
   */
  GeantBasket *GetPBasket() const { return fPBasket.load(); }

  /**
   * @brief Function that set current basket
   * 
   * @param basket Basket that should be set as current
   */
  void SetCBasket(GeantBasket *basket) { fCBasket.store(basket); }

  /**
   * @brief Function that set basket as priority basket 
   * 
   * @param basket Basket that will be stored as priority
   */
  void SetPBasket(GeantBasket *basket) { fPBasket.store(basket); }
#endif

  /**
   * @brief Function that return scheduler
   * @return Scheduler for basket
   */
  GeantScheduler *GetScheduler() const { return fScheduler; }

  /** @brief Function that return name of volume */
  const char *GetName() const { return (fVolume) ? fVolume->GetName() : ClassName(); }

  /**
   * @brief Function that return number assigned to basket
   * @return Number assigned to basket
   */
  Int_t GetNumber() const { return fNumber; }

  /**
   * @brief Function that return volume
   * @return Volume for which applies basket
   */
  TGeoVolume *GetVolume() const { return fVolume; }

  /** @brief Print function */
  virtual void Print(Option_t *option = "") const;

  /**
   * @brief Recycle current Basket
   * 
   * @param b Pointer to current GeantBasket for recycling
   */
  void RecycleBasket(GeantBasket *b);

  /**
   * @brief Function that set feeder queue
   * 
   * @param queue priority_queue for GeantBasket
   */
  void SetFeederQueue(Geant::priority_queue<GeantBasket *> *queue) { fFeeder = queue; }

  /**
   * @brief Function Sizeof() of current basket that had being filled
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

  ClassDef(GeantBasketMgr, 0) // A path in geometry represented by the array of indices
};
#endif
