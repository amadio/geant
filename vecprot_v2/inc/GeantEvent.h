//===--- GeantEvent.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantEvent.h
 * @brief Implementation of event for Geant-V prototype
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_EVENT
#define GEANT_EVENT

#if __cplusplus >= 201103L
#include <atomic>
#endif

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

/** @brief Class GeantEvent that decribes events */
class GeantEvent : public TObject {
private:
#if __cplusplus >= 201103L
  std::atomic_int fEvent;   /** Event number */
  std::atomic_int fSlot;    /** Fixed slot number */
  std::atomic_int fNtracks; /** Number of tracks */
  std::atomic_int fNdone;   /** Number of done tracks */
#else
  // dummy declarations to cheat CINT
  int fEvent;   /** Event number */
  int fSlot;    /** Fixed slot number */
  int fNtracks; /** Number of tracks */
  int fNdone;   /** Number of done tracks */
#endif
  TMutex fMutex; /** mutex for this event */

public:

  /** @brief GeantEvent constructor initialization */
  GeantEvent() : fEvent(0), fSlot(0), fNtracks(0), fNdone(0), fMutex() {}
  
  /**
   * @brief GeantEvent parametrized constructor initialization
   * 
   * @param ntr Number of threads
   */
  GeantEvent(Int_t ntr) : fEvent(0), fSlot(0), fNtracks(ntr), fNdone(0), fMutex() {}
  
  /** @brief GeantEvent destructor */
  ~GeantEvent() {}
  
  /* @brief Function of adding tracks */
  Int_t AddTrack();
#if __cplusplus >= 201103L
  
  /**
   * @brief Function that returns number of event
   * @return Number of event
   */
  Int_t GetEvent() const { return fEvent.load(); }
  
  /**
   * @brief Function that returns fixed number of slot
   * @return Number of slot
   */
  Int_t GetSlot() const { return fSlot.load(); }
  
  /**
   * @brief Function that returns number of tracks
   * @return Number of tracks
   */
  Int_t GetNtracks() const { return fNtracks.load(); }
  
  /**
   * @brief Function that set number of event
   * 
   * @param event Event that should be numbered
   */
  void SetEvent(Int_t event) { fEvent.store(event); }
  
  /**
   * @brief Function that setnumber of slot
   * 
   * @param islot Slot that should be numbered
   */
  void SetSlot(Int_t islot) { fSlot.store(islot); }
  
  /** @brief Reset function */
  void Reset() {
    fNtracks.store(0);
    fNdone.store(0);
  }

  /**
   * @brief Function that check transportation of tracks
   * @return Boolean value that shows if track is transported or not
   */
  Bool_t Transported() const { return ((fNtracks.load() > 0) && (fNtracks == fNdone)); }
  
  /**
   * @brief Function that stop processing of track
   */
  void StopTrack() { fNdone++; }
#else
  
  /**
   * @brief Function that returns number of event
   * @return Number of event
   */
  Int_t GetEvent() const { return fEvent; }
  
  /**
   * @brief Function that returns fixed number of slot
   * @return Number of slot
   */
  Int_t GetSlot() const { return fSlot; }
  
  /**
   * @brief Function that returns number of tracks
   * @return Number of tracks
   */
  Int_t GetNtracks() const { return fNtracks; }
  
  /**
   * @brief Function that set number of event
   * 
   * @param event Event that should be numbered
   */
  void SetEvent(Int_t event) { fEvent = event; }
  
  /**
   * @brief Function that setnumber of slot
   * 
   * @param islot Slot that should be numbered
   */
  void SetSlot(Int_t islot) { fSlot = islot; }
  
  /** @brief Reset function */
  void Reset() { fNtracks = fNdone = 0; }
  
  /**
   * @brief Function that check transportation of tracks
   * @return Boolean value that shows if track is transported or not
   */
  Bool_t Transported() const { return ((fNtracks > 0) && (fNtracks == fNdone)); }
  
 /**
   * @brief Function that stop processing of track
   */
  void StopTrack() {
    fMutex.Lock();
    fNdone++;
    fMutex.UnLock();
  }
#endif

  /** @brief Print function */
  void Print(Option_t *option = "") const;

  ClassDef(GeantEvent, 1) // The G5 event
};
#endif
