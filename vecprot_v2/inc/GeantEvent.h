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

#include <atomic>

/** @brief Class GeantEvent that decribes events */
class GeantEvent {
private:
  bool             fPrioritize;  /** Prioritize this event */
  float            fPriorityThr; /** Priority threshold in percent of max in flight */
  std::atomic_int  fEvent;   /** Event number */
  std::atomic_int  fSlot;    /** Fixed slot number */
  std::atomic_int  fNtracks; /** Number of tracks */
  std::atomic_int  fNdone;   /** Number of done tracks */
  std::atomic_int  fNmax;    /** Maximum number of tracks in flight */
  std::atomic_flag fLock;   /** Lock for priority forcing */
public:

  /** @brief GeantEvent default constructor */
  GeantEvent() : fPrioritize(false), fPriorityThr(0.01), fEvent(0), fSlot(0), fNtracks(0), fNdone(0), fNmax(0), fLock() {}
    
  /** @brief GeantEvent destructor */
  ~GeantEvent() {}
  
  /* @brief Function for accounting adding a new track */
  int AddTrack();
  
  /**
   * @brief Function that returns the event number
   * @return Event number
   */
  int GetEvent() const { return fEvent.load(); }
  
  /**
   * @brief Function that returns the number of slot
   * @return Slot number
   */
  int GetSlot() const { return fSlot.load(); }
  
  /**
   * @brief Function that returns the number of tracks in flight
   * @return Number of tracks in flight
   */
  int GetNinflight() const { return fNtracks.load() - fNdone.load(); }

  /**
   * @brief Function that returns the number of transported tracks
   * @return Number of transported tracks
   */
  int GetNdone() const { return fNdone.load(); }

  /**
   * @brief Function that returns the number of tracks
   * @return Number of tracks
   */
  int GetNtracks() const { return fNtracks.load(); }
  
  /**
   * @brief Function that returns the max number of tracks in flight
   * @return Maximum number of tracks in flight
   */
  int GetNmax() const { return fNmax.load(); }
  
  /**
   * @brief Getter for priority flag
   * @return Priority flag value
   */
  bool IsPrioritized() const { return fPrioritize; }
  
  /**
   * @brief Getter for priority threshold
   * @return Priority flag value
   */
  float GetPriorityThr() const { return fPriorityThr; }
  
  /** @brief Setter for priority threshold */
  void SetPriorityThr(float threshold) { fPriorityThr = threshold; }

  /**
   * @brief Function to set the event number
   * 
   * @param event Event number to be set
   */
  void SetEvent(int event) { fEvent.store(event); }
  
  /**
   * @brief Function to set the slot number
   * 
   * @param islot Slot number to be set
   */
  void SetSlot(int islot) { fSlot.store(islot); }

  /** @brief Prioritize the event */
  bool Prioritize();
 
  /** @brief Reset the event */
  void Reset() {
    fNtracks.store(0);
    fNdone.store(0);
    fNmax.store(0);
    fPrioritize = false;
  }

  /**
   * @brief Function to check if all tracks in the event were transported
   * @return Boolean value that shows if event is fully transported or not
   */
  bool Transported() const { return ((fNtracks.load() > 0) && (fNtracks == fNdone)); }
  
  /**
   * @brief Function to signal that a trach was stopped
   *
   * @return Flag true if stopping qa track started priority mode for the event
   */
  bool StopTrack();

  /** @brief Print function */
  void Print(const char *option = "") const;
};
#endif
