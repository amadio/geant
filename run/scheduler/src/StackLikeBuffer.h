//===--- StackLikeBuffer.h - Geant-V -------------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file StackLikeBuffer.h
 * @brief A track buffer prioritizing high-generation tracks.
 * @details A stack-like buffer holds tracks of different generations in
 *    separate pipelines. New-generation pipelines are flushed into the first
 *    simulation stage earlier than old-generation ones, allowing to transport
 *    with priority newer generations. A special pipeline is reserved for tracks
 *    from events to be prioritized.
 *
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_STACK_LIKE_BUFFER
#define GEANT_STACK_LIKE_BUFFER

#include <VecCore/VecCore>
#include "Geant/Typedefs.h"
#include "Geant/mpmc_bounded_queue.h"
#include "Geant/Track.h"
#include "Geant/Basket.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class TaskData;
class Propagator;
//#include "Geant/Fwd.h"

class StackLikeBuffer {

  using queue_t = mpmc_bounded_queue<Track *>;

protected:
  vector_t<Basket *> fLanes;       ///< Track lanes
  Basket *fPriorityLane = nullptr; ///< Priority lane
  Basket *fStageBuffer  = nullptr; ///< First stage buffer
  int fPriorityEvent    = 0;       ///< Prioritized event
  bool fPriorityMode    = false;   ///< Priority mode
  int fLastLane         = 0;       ///< Last lane containing tracks
  int fNlanes           = 10;      ///< Number of lanes stored
  size_t fNtracks       = 0;       ///< Number of stored tracks

private:
  StackLikeBuffer(const StackLikeBuffer &) = delete;
  StackLikeBuffer &operator=(const StackLikeBuffer &) = delete;

public:
  StackLikeBuffer(int maxgen, TaskData *td) : fNlanes(maxgen + 1)
  {
    if (td->fPropagator->fConfig->fUseNuma) {
      for (int i = 0; i < fNlanes; ++i) {
        fLanes.push_back(new Basket(1000, 0, td->fNode));
      }
      fPriorityLane = new Basket(1000, 0, td->fNode);
    } else {
      for (int i = 0; i < fNlanes; ++i) {
        fLanes.push_back(new Basket(1000, 0));
      }
      fPriorityLane = new Basket(1000, 0);
    }
  }

  ~StackLikeBuffer()
  {
    for (auto basket : fLanes)
      delete basket;
    fLanes.clear();
    delete fPriorityLane;
  }

  /** @brief Add several tracks to the buffer */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void AddTracks(TrackVec_t &tracks)
  {
    for (auto track : tracks)
      AddTrack(track);
  }

  /** @brief Add one track to the buffer */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void AddTrack(Track *track)
  {
    if (fPriorityMode && track->Event() == fPriorityEvent) {
      fPriorityLane->AddTrack(track);
    } else {
      int lane = vecCore::math::Min(track->GetGeneration(), fNlanes - 1);
      fLanes[lane]->AddTrack(track);
      fLastLane = vecCore::math::Max(fLastLane, lane);
    }
    fNtracks++;
  }

  /** @brief Flush a given lane into the stage buffer */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int FlushLane(size_t lane)
  {
    int nflush = fLanes[lane]->size();
    fStageBuffer->AddTracks(fLanes[lane]->Tracks());
    fLanes[lane]->Tracks().clear();
    fNtracks -= nflush;
    return nflush;
  }

  /** @brief Flush the last non-empty lane into the stage buffer */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int FlushLastLane()
  {
    int nflush = FlushLane(fLastLane);
    while (fLastLane > 0 && fLanes[--fLastLane]->size() == 0)
      ;
    return nflush;
  }

  /** @brief Flush the priority lane into the stage buffer */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int FlushPriorityLane()
  {
    int nflush = fPriorityLane->size();
    if (!nflush) return 0;
    fStageBuffer->AddTracks(fPriorityLane->Tracks());
    fPriorityLane->Tracks().clear();
    fNtracks -= nflush;
    return nflush;
  }

  /** @brief Share a number of tracks from lowest generation lanes */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  size_t ShareTracks(size_t ntracks, queue_t &qshare)
  {
    // Start sharing tracks from low lanes
    int lane       = 0;
    size_t nshared = 0;
    while (nshared < ntracks && lane <= fLastLane) {
      size_t nlane = fLanes[lane]->size();
      if (nlane + nshared <= ntracks) {
        // We can share the full lane
        for (auto track : fLanes[lane]->Tracks()) {
          bool shared = qshare.enqueue(track);
          (void)shared;
          assert(shared);
        }
        fLanes[lane]->Tracks().clear();
        fNtracks -= nlane;
        nshared += nlane;
        // continue with next lane
        lane++;
      } else {
        // We can share only ntracks - nshared tracks
        while (nshared < ntracks) {
          auto track = fLanes[lane]->Tracks().back();
          fLanes[lane]->Tracks().pop_back();
          bool shared = qshare.enqueue(track);
          (void)shared;
          assert(shared);
          fNtracks--;
          nshared++;
        }
        return nshared;
      }
    }
    return nshared;
  }

  /** @brief Setter for the first stage buffer */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void SetStageBuffer(Basket *buffer) { fStageBuffer = buffer; }

  /** @brief Getter for the priority event */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int GetPriorityEvent() const { return fPriorityEvent; }

  /** @brief Setter for the priority event */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void SetPriorityEvent(int event) { fPriorityEvent = event; }

  /** @brief Getter for the priority mode */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  bool IsPrioritized() const { return fPriorityMode; }

  /** @brief Setter for the priority mode */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void SetPriorityMode(bool flag) { fPriorityMode = flag; }

  /** @brief Getter for number of stacked tracks */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int GetNtracks() const { return fNtracks; }

  /** @brief Getter for number of stacked tracks */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int GetNprioritized() const { return fPriorityLane->size(); }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
