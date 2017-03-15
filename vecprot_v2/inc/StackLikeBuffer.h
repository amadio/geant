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
#include "GeantTrack.h"
#include "Basket.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class GeantPropagator;
//#include "GeantFwd.h"

class StackLikeBuffer {

protected:  
  vector_t<Basket *> fLanes;       ///< Track lanes
  Basket *fPriorityLane = nullptr; ///< Priority lane
  Basket *fStageBuffer = nullptr;  ///< First stage buffer
  int fPriorityEvent = 0;          ///< Prioritized event
  bool fPriorityMode = false;      ///< Priority mode
  int fLastLane = 0;               ///< Last lane containing tracks
  int fNlanes = 10;                ///< Number of lanes stored
  
private:
  StackLikeBuffer(const StackLikeBuffer &) = delete;
  StackLikeBuffer &operator=(const StackLikeBuffer &) = delete;

public:
  StackLikeBuffer(int maxgen, GeantTaskData *td) : fNlanes(maxgen+1) 
  {
    for (int i=0; i<fNlanes; ++i) {
      fLanes.push_back(new Basket(1000, 0, td->fPropagator->fNuma));
    }
    fPriorityLane = new Basket(1000, 0, td->fPropagator->fNuma);
  }
  
  ~StackLikeBuffer()
  {
     for (auto basket : fLanes) delete basket;
     fLanes.clear();
     delete fPriorityLane;
  }
  
  /** @brief Add several tracks to the buffer */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void AddTracks(TrackVec_t &tracks) { 
    for (auto track : tracks) 
      AddTrack(track);
  }

  /** @brief Add one track to the buffer */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void AddTrack(GeantTrack *track)
  {
    if (fPriorityMode && track->fEvent == fPriorityEvent) {
      fPriorityLane->AddTrack(track);
    } else {
      int lane = vecCore::math::Min(track->fGeneration, fNlanes-1);
      fLanes[lane]->AddTrack(track);
      fLastLane = vecCore::math::Max(fLastLane, lane);
    }
  }
  
  /** @brief Flush a given lane into the stage buffer */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int FlushLane(size_t lane)
  { 
    int nflush = fLanes[lane]->size();
    fStageBuffer->AddTracks(fLanes[lane]->Tracks());
    fLanes[lane]->Tracks().clear();
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
  int FlushPriorityLane() {
    int nflush = fPriorityLane->size();
    if (!nflush) return 0;
    fStageBuffer->AddTracks(fPriorityLane->Tracks());
    fPriorityLane->Tracks().clear();
    return nflush;
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
  int GetNtracks() const 
  { 
    int ntracks = 0;
    for (int lane = fLastLane; lane >= 0; --lane)
      ntracks += fLanes[lane]->size();
    return ntracks;
  }

  /** @brief Getter for number of stacked tracks */
  GEANT_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  int GetNprioritized() const { return fPriorityLane->size(); }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
