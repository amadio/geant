//===--- GeantEvent.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantEvent.h
 * @brief Implementation of event for GeantV prototype
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
#endif
  TMutex fMutex; /** mutex for this event */

public:

  /** @brief GeantEvent default constructor */
  GeantEvent() : fEvent(0), fSlot(0), fNtracks(0), fNdone(0), fMutex() {}
    
  /** @brief GeantEvent destructor */
  ~GeantEvent() {}
  
  /* @brief Function for accounting adding a new track */
  Int_t AddTrack();
#if __cplusplus >= 201103L
  
  /**
   * @brief Function that returns the event number
   * @return Event number
   */
  Int_t GetEvent() const { return fEvent.load(); }
  
  /**
   * @brief Function that returns the number of slot
   * @return Slot number
   */
  Int_t GetSlot() const { return fSlot.load(); }
  
  /**
   * @brief Function that returns the number of transported tracks
   * @return Number of transported tracks
   */
  Int_t GetNdone() const { return fNdone.load(); }

  /**
   * @brief Function that returns the number of tracks
   * @return Number of tracks
   */
  Int_t GetNtracks() const { return fNtracks.load(); }
  
  /**
   * @brief Function to set the event number
   * 
   * @param event Event number to be set
   */
  void SetEvent(Int_t event) { fEvent.store(event); }
  
  /**
   * @brief Function to set the slot number
   * 
   * @param islot Slot number to be set
   */
  void SetSlot(Int_t islot) { fSlot.store(islot); }
  
  /** @brief Reset the event */
  void Reset() {
    fNtracks.store(0);
    fNdone.store(0);
  }

  /**
   * @brief Function to check if all tracks in the event were transported
   * @return Boolean value that shows if event is fully transported or not
   */
  Bool_t Transported() const { return ((fNtracks.load() > 0) && (fNtracks == fNdone)); }
  
  /**
   * @brief Function to signal that a trach was stopped
   */
  void StopTrack() { fNdone++; }
#endif

  /** @brief Print function */
  void Print(Option_t *option = "") const;

  ClassDef(GeantEvent, 1) // The G5 event
};
#endif
