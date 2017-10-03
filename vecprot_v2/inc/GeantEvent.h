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
#include <vector>
#include "Geant/Config.h"
#include "base/Vector3D.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTrack;
class GeantRunManager;
class GeantTaskData;

/** @brief Class GeantEvent that decribes events */
class GeantEvent {

private:
  vecgeom::Vector3D<double> fVertex;     /** Vertex position */
  bool             fPrioritize = false;  /** Prioritize this event */
  bool             fTransported = false; /** Event transported */
  float            fPriorityThr = 0.01; /** Priority threshold in percent of max in flight */
  int              fEvent = 0;      /** Event number */
  int              fSlot = 0;       /** Fixed slot number === to be removed ===*/
  int              fNprimaries = 0; /** Number of primaries */
  std::atomic_int  fNtracks;        /** Number of tracks */
  std::atomic_int  fNdone;          /** Number of done tracks */
  std::atomic_int  fNmax;           /** Maximum number of tracks in flight */
  std::atomic_flag fLock;           /** Lock for priority forcing */
  std::vector<GeantTrack*> fPrimaries; /** Vector containing all primary tracks */
public:

  std::atomic_int  fNfilled;        /** Number of tracks copied in buffer */
  std::atomic_int  fNdispatched;    /** Number of tracks dispatched */

  /** @brief GeantEvent default constructor */
  GeantEvent() : fNtracks(0), fNdone(0), fNmax(0), fLock(), fNfilled(0), fNdispatched(0) {}

  /** @brief GeantEvent destructor */
  ~GeantEvent() {}

  /* @brief Function for accounting adding a new track */
  int AddTrack();

  /* @brief Function for accounting adding a new track */
  int AddPrimary(GeantTrack *track) { fPrimaries.push_back(track); return AddTrack(); }

  /* @brief Crear the event and release all primaries */
  void Clear();

  /** @brief Dispatch track. */
  GEANT_FORCE_INLINE
  int DispatchTrack(bool &valid) {
    int itr = fNdispatched.fetch_add(1);
    valid = itr < fNprimaries;
    return itr;
  }

  /** @brief Check if event is dispatched. */
  GEANT_FORCE_INLINE
  bool IsDispatched() const { return (fNdispatched.load() >= fNprimaries); }

  /* @brief Function for retrieving a primary. No range check. */
  GEANT_FORCE_INLINE
  GeantTrack *GetPrimary(int i) { return fPrimaries[i]; }

  /* @brief Function for retrieving a primary. No range check. */
  GEANT_FORCE_INLINE
  int GetNprimaries() const { return fNprimaries; }

  /* @brief Function for retrieving a primary. No range check. */
  GEANT_FORCE_INLINE
  void SetNprimaries(int nprim) { fNprimaries = nprim; fPrimaries.reserve(nprim);}

  /** @brief Function that returns the event vertex */
  GEANT_FORCE_INLINE
  vecgeom::Vector3D<double> GetVertex() const { return fVertex; }

  /** @brief Function that returns the event number */
  GEANT_FORCE_INLINE
  int GetEvent() const { return fEvent; }

  /** @brief Function that returns the slot number */
  GEANT_FORCE_INLINE
  int GetSlot() const { return fSlot; }

  /**
   * @brief Function that returns the number of tracks in flight
   * @return Number of tracks in flight
   */
  GEANT_FORCE_INLINE
  int GetNinflight() const { return fNtracks.load() - fNdone.load(); }

  /**
   * @brief Function that returns the number of transported tracks
   * @return Number of transported tracks
   */
  GEANT_FORCE_INLINE
  int GetNdone() const { return fNdone.load(); }

  /**
   * @brief Function that returns the number of tracks
   * @return Number of tracks
   */
  GEANT_FORCE_INLINE
  int GetNtracks() const { return fNtracks.load(); }

  /**
   * @brief Function that returns the max number of tracks in flight
   * @return Maximum number of tracks in flight
   */
  GEANT_FORCE_INLINE
  int GetNmax() const { return fNmax.load(); }

  /**
   * @brief Getter for priority flag
   * @return Priority flag value
   */
  GEANT_FORCE_INLINE
  bool IsPrioritized() const { return fPrioritize; }

  /**
   * @brief Getter for priority threshold
   * @return Priority flag value
   */
  GEANT_FORCE_INLINE
  float GetPriorityThr() const { return fPriorityThr; }

  /** @brief Setter for priority threshold */
  GEANT_FORCE_INLINE
  void SetPriorityThr(float threshold) { fPriorityThr = threshold; }

  /**
   * @brief Function to set the event number
   *
   * @param event Event number to be set
   */
  GEANT_FORCE_INLINE
  void SetEvent(int event) { fEvent = event; }

  /** @brief Function to set the vertex */
  GEANT_FORCE_INLINE
  void SetVertex(double x, double y, double z) { fVertex.Set(x, y, z); }

  /**
   * @brief Function to set the slot number
   *
   * @param islot Slot number to be set
   */
  GEANT_FORCE_INLINE
  void SetSlot(int islot) { fSlot = islot; }

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
  GEANT_FORCE_INLINE
  bool Transported() const { return fTransported; }

  /**
   * @brief Function to signal that a trach was stopped
   *
   * @return Flag true if stopping qa track started priority mode for the event
   */
  bool StopTrack(GeantRunManager *runmgr, GeantTaskData *td);

  /** @brief Print function */
  void Print(const char *option = "") const;
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
