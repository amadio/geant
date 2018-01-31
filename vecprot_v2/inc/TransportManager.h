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

#include <algorithm>
#include "Geant/Typedefs.h"

#include "Geant/Config.h"
#include "Geant/math_wrappers.h"
#include "GeantTrack.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class Basket;

/**
 * @brief Class GeantTrack
 */
namespace TransportManager {

  /**
   * @brief Check if the geometry location changed for a vector of tracks
   *
   * @param tracks Vector of tracks (AOS)
   * @param ntracks Number of tracks
   * @param td TaskData object
   */
  VECCORE_ATT_HOST_DEVICE
  int CheckSameLocation(TrackVec_t &tracks,
                         int ntracks,
                         GeantTaskData *td);

  /**
   * @brief Check if the geometry location changed for a track
   *
   * @param track Track reference
   * @param td TaskData object
   */
  VECCORE_ATT_HOST_DEVICE
  int CheckSameLocationSingle(GeantTrack &track,
                         GeantTaskData *td);

  /**
   * @brief Compute transport length for a vector of tracks
   *
   * @param tracks Vector of tracks (AOS)
   * @param ntracks Number of tracks
   * @param td TaskData object
   */
  VECCORE_ATT_HOST_DEVICE
  void ComputeTransportLength(TrackVec_t &tracks,
                              int ntracks,
                              GeantTaskData *td);

  /**
   * @brief Compute transport length for single track
   *
   * @param track Track reference
   * @param td TaskData object
   */
  VECCORE_ATT_HOST_DEVICE
  void ComputeTransportLengthSingle(GeantTrack &track,
                                    GeantTaskData *td);

  /**
   * @brief Function that provides postponed action for tracks
   *
   * @param ntracks Number of tracks
   */  
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  TransportAction_t PostponedAction(int ntracks) {
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
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void MoveTrack(int itr, TrackVec_t &input, TrackVec_t &output) {
#ifndef VECCORE_CUDA
    auto it = input.begin() + itr;
    std::move(it, it+1, std::back_inserter(output));
#else
    output.push_back(input[itr]);
#endif
    input.erase(input.begin() + itr);
  }

  /**
   * @brief Function rotating a track to the end of the vector
   *
   * @param itr track index
   * @param vec Vector of input tracks
   */  
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void RotateTrack(int itr, TrackVec_t &vec) {
#ifndef VECCORE_CUDA
    std::rotate(vec.begin()+itr, vec.begin()+itr+1, vec.end());
#else
    auto temp = vec[itr];
    auto it = vec.begin() + itr;
    vec.erase(it);
    vec.push_back(temp);
#endif
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
  void PropagateInVolumeSingle(GeantTrack &track,
                               double crtstep,
                               GeantTaskData *td);

  /**
   * @brief Propagate a vector of tracks according their proposed physics steps
   *
   * @param tracks Vector of tracks (AOS) to be propagated
   * @param td Task data object
   */
  VECCORE_ATT_HOST_DEVICE
  int PropagateTracks(TrackVec_t &tracks,
                      GeantTaskData *td);

  /**
   * @brief Propagate a vector of tracks in scalar mode
   *
   * @param tracks Container of tracks
   * @param td Task data object
   * @param stage Transport is done in several stages, can resume from a given one
   */
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
  int PropagateSingleTrack(TrackVec_t &tracks,
                           int &itr,
                           GeantTaskData *td,
                           int stage);

  VECCORE_ATT_HOST_DEVICE
  int PropagateSingleTrack(GeantTrack *track,
                           Basket *output,
                           GeantTaskData *td,
                           int stage);

  /** @brief Function that returns safe length */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double SafeLength(const GeantTrack &track, double Bz, double eps = 1.E-4) {
    // Returns the propagation length in field such that the propagated point is
    // shifted less than eps with respect to the linear propagation.
    return 2. * sqrt(eps / track.Curvature(Bz));
  }
  
  /** @brief Function allowing to debug a step */
  VECCORE_ATT_HOST_DEVICE
 bool BreakOnStep(TrackVec_t &tracks, int evt, int trk, int stp, int nsteps, const char *msg, int itr=-1);

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
