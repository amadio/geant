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
  
private:
  size_t fBlockSize;      // Number of tracks stored by each block
  int fNode;              // Numa id for the managed blocks
  NumaBlockMgr<GeantTrack> fBlockMgr; // Block manager for tracks
    
public:
  /** @brief Constructor */
  TrackManager(size_t nblocks, size_t block_size, int numa_node=0)
    :fBlockSize(block_size), fNode(numa_node), fBlockMgr(numa_node, block_size) {
  
    for (size_t i=0; i<nblocks; ++i) {
      NumaBlock<GeantTrack> *block = fBlockMgr.AddBlockAndRegister();
      for (size_t itr=0; itr<fBlockSize; ++itr) block[itr].fBindex = itr;
    }
  }

  /** @brief Destructor */
  ~TrackManager() = delete;
  
  /** @brief Returns a reference of a track from the container */
  GEANT_INLINE
  GeantTrack & const GetTrack() { return fBlockMgr.GetObject(); }
  
  /** @brief Service to get the NUMA block a track belongs to */
  GEANT_INLINE
  static
  NumaBlock<GeantTrack> *GetBlock(GeantTrack & const track) {
    // The magic here is that tracks are allocated contiguous in blocks and the
    // block address is stored just before the array of tracks. So we need to
    // store just the track index in the track structure to retrieve the address
    // of the block
    return ( *((NumaBlock<GeantTrack>**)(track - track.fBindex) - 1) );
  }

  /** @brief Release a track from its block */
  GEANT_INLINE
  bool ReleaseTrack(GeantTrack & const track) {
    return ( fBlockMgr.ReleaseObject(GetBlock(track)) );
  }
  
  /** @brief Check if track belongs to container */
  GEANT_INLINE
  bool OwnsTrack(GeantTrack *track) const {
    return ( (size_t)(track - fTracks) < fSize*sizeof(GeantTrack) );
  } 

};
} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
