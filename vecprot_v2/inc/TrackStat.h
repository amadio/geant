//===--- TrackStat.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TrackStat.h
 * @brief Statistics class for tracks.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TRACK_STAT
#define GEANT_TRACK_STAT

#include "Geant/Typedefs.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;

class TrackStat {

private:
  GeantTaskData *fTd = nullptr;   ///< Task data
  int fNtotal = 0;                ///< Total number of tracks injected by the thread
  int fNstacked = 0;              ///< Number of tracks in the stack-like buffer
  int fNstagebuff = 0;            ///< Number of tracks in stage buffers
  int fNbasketized = 0;           ///< Number of tracks in un-flushed baskets

public:
  /** @brief Track stat constructor */
  VECCORE_ATT_HOST_DEVICE
  TrackStat(GeantTaskData *td) : fTd(td) {}

  /** @brief Track stat destructor */
  VECCORE_ATT_HOST_DEVICE
  ~TrackStat() {}

  /** @brief Clear stats */
  VECCORE_ATT_HOST_DEVICE
  void AddTracks(int ntracks) { fNtotal += ntracks; }

  /** @brief Clear stats */
  VECCORE_ATT_HOST_DEVICE
  void Clear() { fNtotal = 0; fNstacked = 0; fNstagebuff = 0; fNbasketized = 0; }
  
  /** @brief Make the count balance for tracks */
  int CountBalance();
  
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
