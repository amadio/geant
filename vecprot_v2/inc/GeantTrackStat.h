//===--- GeantTrackStat.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTrackStat.h
 * @brief Implementation of statistics for track array in Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TRACK_STAT
#define GEANT_TRACK_STAT

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

class GeantTrack;
class GeantTrack_v;

/**
 * @brief Class GeantTrackStat
 * @details Used statistic object for a track array
 * 
 */
class GeantTrackStat : public TObject {
public:
  Int_t fNslots;   /** Number of event slots */
  Int_t *fNtracks; /** [fNslots] Number of tracks from an event */
  Int_t *fNsteps;  /**[fNslots] Cumulated number of steps per event */
  TMutex fMutex;   /** Mutex */

private:

  /**
   * @brief Copy constructor
   */
  GeantTrackStat(const GeantTrackStat &other);

  /**
   * @brief Operator =
   * 
   * @param other ?????
   */
  GeantTrackStat &operator=(const GeantTrackStat &other);

public:

  /** @brief GeantTrackStat constructor */
  GeantTrackStat() : TObject(), fNslots(0), fNtracks(0), fNsteps(0), fMutex() {}

  /**
   * @brief GeantTrackStat constructor
   * 
   * @param nslots Number of event slots
   */
  GeantTrackStat(Int_t nslots);

  /** @brief GeantTrackStat destructor */
  virtual ~GeantTrackStat();

  /**
   * @brief Operator +=
   * 
   * @param other ??????
   */
  GeantTrackStat &operator+=(const GeantTrackStat &other);

  /**
   * @brief Operator -=
   * 
   * @param other ?????
   */
  GeantTrackStat &operator-=(const GeantTrackStat &other);

  /**
   * @brief Function that add track
   * 
   * @param track Track that should be added
   */
  void AddTrack(const GeantTrack &track);

  /**
   * @brief Function that add track from track_v array
   * 
   * @param trackv Track that should be added
   * @param itr Track ID ?????
   */
  void AddTrack(const GeantTrack_v &trackv, Int_t itr);

  /**
   * @brief Function that add tracks from track_v array
   * 
   * @param trackv Tracks that should be added
   */
  void AddTracks(const GeantTrack_v &trackv);

  /**
   * @brief Function that remove tracks
   * 
   * @param trackv Tracks that should be deleted
   */
  void RemoveTracks(const GeantTrack_v &trackv);

  /**
   * @brief Function that init array of event slots ?????
   * 
   * @param nslots Number of event slots
   */
  void InitArrays(Int_t nslots);

  /** @brief Print function */
  void Print(Option_t *option = "") const;

  /** @brief Reset function */
  void Reset();

  ClassDef(GeantTrackStat, 1) // Track statistics
};
#endif
