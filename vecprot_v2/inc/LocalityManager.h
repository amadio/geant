//===--- LocalityManager.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file LocalityManager.h
 * @brief Manager class for handling locality.
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_LOCALITY_MANAGER
#define GEANT_LOCALITY_MANAGER

#include "Geant/Config.h"
#include "TrackManager.h"
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
class LocalityManager {
  using size_t = std::size_t;
  
private:
  static LocalityManager *fgInstance; /** Singleton instance */
  bool fInitialized;               /** Flag for initialization */
  NumaPolicy fPolicy;              /** Numa policy to be used */
  int fNnodes;                     /** Number of nodes */
  TrackManager **fTrackMgr;        /** [fNnodes] Array of track managers */
  size_t fNblocks;                 /** Number of initial track blocks */
  size_t fBlockSize;               /** Number of tracks stored by each block */
  int    fMaxdepth;                /** Maximum geometry depth */

private:
  /** @brief Constructor */
  LocalityManager();    
public:
 /** @brief Function that creates locality manager instance **/
  static
  LocalityManager* Instance();

  /** @brief Destructor */
  ~LocalityManager();
  
  /** @brief Setter for the number of locality nodes.*/
  GEANT_FORCE_INLINE
  void SetNnodes(int nnodes) { fNnodes = nnodes; }

  /** @brief Getter for the number of locality nodes.*/
  GEANT_FORCE_INLINE
  int GetNnodes() const { return fNnodes; }

  /** @brief Getter for the number of queued blocks.*/
  GEANT_FORCE_INLINE
  int GetNqueued() const
  {
    int nqueued = 0;
    for (auto i=0; i<fNnodes; ++i) nqueued += fTrackMgr[i]->GetNqueued();
    return nqueued;
  }

  /** @brief Getter for the number of queued blocks.*/
  GEANT_FORCE_INLINE
  int GetNallocated() const
  {
    int nblocks = 0;
    for (auto i=0; i<fNnodes; ++i) nblocks += fTrackMgr[i]->GetNblocks();
    return nblocks;
  }

  /** @brief Getter for the number of released blocks.*/
  GEANT_FORCE_INLINE
  int GetNreleased() const
  {
    int nblocks = 0;
    for (auto i=0; i<fNnodes; ++i) nblocks += fTrackMgr[i]->GetNreleased();
    return nblocks;
  }

  /** @brief Setter for the locality policy.*/
  GEANT_FORCE_INLINE
  void SetPolicy(NumaPolicy::EPolicyType policy) { fPolicy.SetPolicy(policy); }

  /** @brief Getter for the locality policy.*/
  GEANT_FORCE_INLINE
  NumaPolicy &GetPolicy() const { return (NumaPolicy&)fPolicy; }
  
  /** @brief Setter for the number of blocks to allocate.*/
  GEANT_FORCE_INLINE
  void SetNblocks(size_t nblocks) { fNblocks = nblocks; }

  /** @brief Getter for the number of blocks to allocate.*/
  GEANT_FORCE_INLINE
  int GetNblocks() const { return fNblocks; }
  
  /** @brief Setter for the block size.*/
  GEANT_FORCE_INLINE
  void SetBlockSize(size_t size) { fBlockSize = size; }

  /** @brief Getter for the block size.*/
  GEANT_FORCE_INLINE
  int GetBlockSize() const { return fBlockSize; }

  /** @brief Setter for the maximum geometry depth.*/
  GEANT_FORCE_INLINE
  void SetMaxDepth(int maxdepth) { fMaxdepth = maxdepth; }

  /** @brief Getter for the maximum geometry depth.*/
  GEANT_FORCE_INLINE
  int GetMaxDepth() const { return fMaxdepth; }
  
  /** @brief Initialize locality manager and allocate data.*/
  void Init();
  
  /** @brief Getter for the initialization flag.*/
  GEANT_FORCE_INLINE
  bool IsInitialized() const { return fInitialized; }

  /** @brief Getter for track managers per locality node.*/
  GEANT_FORCE_INLINE
  TrackManager &GetTrackManager(int inode) const { return *fTrackMgr[inode]; }
      
  /** @brief Service to recycle tracks */
  GEANT_FORCE_INLINE
  bool ReleaseTrack(GeantTrack const &track) {
    int node = fTrackMgr[0]->GetNode(track);
    return ( fTrackMgr[node]->ReleaseTrack(track) );
  }

};
} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
