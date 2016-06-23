//===--- TransportManager.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TransportManager.h
 * @brief Namespace providing methods for querying and stepping in geometry
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TRANSPORT_MANAGER
#define GEANT_TRANSPORT_MANAGER

#include <vector>
#include <algorithm>

#include "Geant/Config.h"
#include "Geant/Math.h"
#include "GeantTrack.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;

/**
 * @brief Class GeantTrack
 */
namespace TransportManager {
  typedef std::vector<GeantTrack * const> TrackVec_t;

  /**
   * @brief Check if the geometry location changed for a vector of tracks
   *
   * @param tracks Vector of tracks (AOS)
   * @param ntracks Number of tracks
   * @param td TaskData object
   */
  GEANT_CUDA_BOTH_CODE
  int CheckSameLocation(TrackVec_t &tracks,
                         int ntracks,
                         GeantTaskData *td);

  /**
   * @brief Check if the geometry location changed for a track
   *
   * @param track Track reference
   * @param td TaskData object
   */
  GEANT_CUDA_BOTH_CODE
  int CheckSameLocationSingle(GeantTrack &track,
                         GeantTaskData *td);

  /**
   * @brief Compute transport length for a vector of tracks
   *
   * @param tracks Vector of tracks (AOS)
   * @param ntracks Number of tracks
   * @param td TaskData object
   */
  GEANT_CUDA_BOTH_CODE
  void ComputeTransportLength(TrackVec_t &tracks,
                              int ntracks,
                              GeantTaskData *td);

  /**
   * @brief Compute transport length for single track
   *
   * @param track Track reference
   * @param td TaskData object
   */
  GEANT_CUDA_BOTH_CODE
  void ComputeTransportLengthSingle(GeantTrack &track,
                                    GeantTaskData *td);

  /**
   * @brief Function that provides postponed action for tracks
   *
   * @param ntracks Number of tracks
   */  
  GEANT_CUDA_BOTH_CODE
  GEANT_INLINE
  TransportAction_t PostponedAction(int ntracks) const {
    // Check the action to be taken according the current policy
    const int kMinVecSize = 4; // this should be retrieved from elsewhere
    if (!ntracks)
      return kDone;
    if (ntracks < kMinVecSize)
      return kSingle; // kPostpone can be an alternative policy
    return kVector;
  }

  /**
   * @brief Function moving a track from one vector to the end of another
   *
   * @param itr track index
   * @param input Vector of input tracks
   * @param output Receiver vector
   */  
  GEANT_CUDA_BOTH_CODE
  GEANT_INLINE
  void MoveTrack(int itr, TrackVec_t &input, TrackVec_t &output) {
    std::move(input.begin()+itr, input.begin()+itr+1, std::back_inserter(output));
    input.erase(input.begin() + itr);
  }

  /**
   * @brief Function rotating a track to the end of the vector
   *
   * @param itr track index
   * @param vec Vector of input tracks
   */  
  GEANT_CUDA_BOTH_CODE
  GEANT_INLINE
  void RotateTrack(int itr, TrackVec_t &vec) {    
    std::rotate(vec.begin()+itr, vec.begin()+itr+1, vec.end());
  }
  
  /**
   * @brief Function that postpones propagation of tracks
   *
   * @param input Tracks to be postponed
   * @param output New vector containing postponed tracks
   */
  int PostponeTracks(TrackVec_t &input, TrackVec_t &output);
  
  /**
   * @brief Propagate a vector of tracks in their current volume
   *
   * @param tracks Vector of tracks (AOS) to be propagated
   * @param ntracks Number of tracks
   * @param crtstep Array of steps to propagate with
   * @param td Task data object
   */
  GEANT_CUDA_BOTH_CODE
  void PropagateInVolume(TrackVec_t &tracks,
                         int ntracks,
                         const double *crtstep,
                         GeantTaskData *td);

  /**
   * @brief Propagate a single track in its current volume
   *
   * @param track Track to be propagated
   * @param crtstep Step to propagate with
   * @param td Task data object
   */
  GEANT_CUDA_BOTH_CODE
  void PropagateInVolumeSingle(GeantTrack &track,
                               double crtstep,
                               GeantTaskData *td);

  /**
   * @brief Propagate a vector of tracks according their proposed physics steps
   *
   * @param tracks Vector of tracks (AOS) to be propagated
   * @param td Task data object
   */
  GEANT_CUDA_BOTH_CODE
  int PropagateTracks(TrackVec_t &tracks,
                      GeantTaskData *td);

  /**
   * @brief Propagate a vector of tracks in scalar mode
   *
   * @param tracks Container of tracks
   * @param td Task data object
   * @param stage Transport is done in several stages, can resume from a given one
   */
  GEANT_CUDA_BOTH_CODE
  int PropagateTracksScalar(TrackVec_t &tracks,
                            GeantTaskData *td,
                            int stage = 0);

  /**
   * @brief Propagate a single track according to the proposed physics step
   *
   * @param tracks Container of tracks
   * @param itr Track number
   * @param td Task data object
   * @param stage Transport is done in several stages, can resume from a given one
   */
  GEANT_CUDA_BOTH_CODE
  int PropagateSingleTrack(TrackVec_t &tracks,
                           int &itr,
                           GeantTaskData *td,
                           int stage);

  /** @brief Function that returns safe length */
  GEANT_CUDA_BOTH_CODE
  GEANT_INLINE
  double SafeLength(const GeantTrack &track, double eps = 1.E-4) {
    // Returns the propagation length in field such that the propagated point is
    // shifted less than eps with respect to the linear propagation.
    return 2. * sqrt(eps / track.Curvature());
  }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
