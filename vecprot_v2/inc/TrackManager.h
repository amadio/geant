//===--- TrackManager.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TrackManager.h
 * @brief Track manager for a given NUMA node
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_TRACK_MANAGER
#define GEANT_TRACK_MANAGER

#include "Geant/Config.h"
#include "NumaBlockMgr.h"
#include "GeantTrack.h"

/**
 * @brief Class TrackManager
 * @detailed The class is managing allocation of blocks of tracks for a given NUMA
 *           node. The allocation for single blocks is using numa_aligned_malloc.
 */
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
/**
 * @brief Class TrackManager
 * @detailed The track manager handles allocation and cleaning of track blocks.
 */
class TrackManager {
  using size_t = std::size_t;
  using NumaTrackBlock_t = NumaBlock<GeantTrack>;
  
private:
//  size_t fBlockSize;      // Number of tracks stored by each block
//  int fNode;              // Numa id for the managed blocks
  NumaBlockMgr<GeantTrack> fBlockMgr; // Block manager for tracks
    
public:
  /** @brief Constructor */
  TrackManager(size_t nblocks, size_t block_size, int numa_node=0)
              : fBlockMgr(nblocks, numa_node, block_size) {      
//    std::cout << "Track manager for numa #" << numa_node << " created at:" << this << std::endl;
//    std::cout << "   NumaBlockMgr: " << &fBlockMgr << std::endl;
  }

  /** @brief Destructor */
  ~TrackManager() {}
  
  /** @brief Returns a reference of a track from the container */
  GEANT_FORCE_INLINE
  GeantTrack &GetTrack() { 
    size_t index;
    GeantTrack &track = fBlockMgr.GetObject(index);
    //track.Clear();
    track.fBindex = index;
    return track;
  }

  /** @brief Get number of queued blocks */
  int GetNqueued() const { return fBlockMgr.GetNqueued(); }

  /** @brief Get number of queued blocks */
  int GetNblocks() const { return fBlockMgr.GetNblocks(); }

  /** @brief Get number of queued blocks */
  int GetNreleased() const { return fBlockMgr.GetNreleased(); }

  /** @brief Service to get the NUMA block a track belongs to */
  GEANT_FORCE_INLINE
  static
  NumaTrackBlock_t *GetBlock(GeantTrack const &track) {
    // The magic here is that tracks are allocated contiguous in blocks and the
    // block address is stored just before the array of tracks. So we need to
    // store just the track index in the track structure to retrieve the address
    // of the block
    auto tr_size = GeantTrack::SizeOfInstance();
    return *((NumaTrackBlock_t**)((char*)&track - tr_size*track.fBindex - sizeof(NumaTrackBlock_t*) - sizeof(GeantTrack*)));
  }

  /** @brief Get a new NUMA block */
  GEANT_FORCE_INLINE
  NumaTrackBlock_t *GetNewBlock() {
    return ( fBlockMgr.GetBlock() );
  }

  /** @brief Recycle NUMA block */
  GEANT_FORCE_INLINE
  void RecycleBlock(NumaTrackBlock_t *block) {
    fBlockMgr.RecycleBlock(block);
  }
 
  /** @brief Service to get the NUMA node corresponding to a track */
  GEANT_FORCE_INLINE
  int GetNode(GeantTrack const &track) {
    return ( TrackManager::GetBlock(track)->GetNode() );
  }

  /** @brief Release a track from its block */
  GEANT_FORCE_INLINE
  bool ReleaseTrack(GeantTrack const &track) {
    NumaTrackBlock_t *block = GetBlock(track);
//    std::cout << " ... releasing track " << &track << " from block " << block << std::endl;
    return ( fBlockMgr.ReleaseObject(block) );
  }
  
};
} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
