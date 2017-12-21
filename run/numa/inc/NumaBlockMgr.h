//===--- NumaBlockMgr.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaBlockMgr.h
 * @brief Memory manager for NUMA-aware data blocks
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_BLOCK_ALLOCATOR
#define GEANT_NUMA_BLOCK_ALLOCATOR

#include <atomic>
#include <vector>
#include <cassert>
#include <iostream>
#include <type_traits>
#include "Geant/Config.h"
#include "mpmc_bounded_queue.h"
#include "NumaBlock.h"

/**
 * @brief Class NumaBlockMgr
 * @detailed The class is managing allocation of blocks of POD data for a given 
 *           NUMA node. The allocation for single blocks is using numa_aligned_malloc.
 */
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

template <typename T> class NumaBlockMgr {
  using size_t = std::size_t;
  using numa_block_ptr = NumaBlock<T>*;
  using queue_t = mpmc_bounded_queue<numa_block_ptr>;

  static size_t const cacheline_size = 64;
  static size_t const queue_buff_size = 1<<16;
  typedef char cacheline_pad_t[cacheline_size];

private:
  std::atomic<numa_block_ptr> fCurrent; // Current block being distributed

  cacheline_pad_t pad0_;   //! Padding to protect the other data from the hot cache line above
  
  int             fNode;  // Numa node id
  std::atomic_int fNblocks; // Number of blocks in flight
  std::atomic_int fNreleased; // Number of blocks in flight
  size_t          fBlockSize;  // Numa block size
  queue_t         fBlocks;  // Queue of free blocks

public:

  /** @brief Constructor providing number of blocks to be initially created */
  NumaBlockMgr(size_t nblocks, int numa_node, size_t bsize) 
    : fCurrent(nullptr), fNode(numa_node), fNblocks(0), fNreleased(0),
      fBlockSize(bsize), fBlocks(queue_buff_size)
  {
    // Constructor creating nblocks
    numa_block_ptr block;
    for (size_t i=0; i<nblocks; ++i) {
      int id = fNblocks.load();
      fNblocks++;
      block = NumaBlock<T>::MakeInstance(fBlockSize, fNode, id);
      if (i == 0)
        fCurrent.store(block);
      else
        fBlocks.enqueue(block);
    }
  }
  
  /** @brief Add a free block */
  numa_block_ptr AddBlock()
  {
    int id = fNblocks.load();
    fNblocks++;
    return ( NumaBlock<T>::MakeInstance(fBlockSize, fNode, id) );
  }

  /** @brief Get number of queued blocks */
  int GetNqueued() const { return fBlocks.size(); }

  /** @brief Get number of blocks in flight */
  int GetNblocks() const { return fNblocks.load(); }

  /** @brief Get number of released */
  int GetNreleased() const { return fNreleased.load(); }

  /** @brief Destructor*/
  ~NumaBlockMgr() {
    numa_block_ptr block;
//    std::cout << "deleting block manager: " << this << std::endl;
    while (fBlocks.dequeue(block)) {
      if (block == fCurrent.load()) continue;
      NumaBlock<T>::ReleaseInstance(block);
    }
//    std::cout << "  delete fCurrent\n";
    NumaBlock<T>::ReleaseInstance(fCurrent.load());
  }

  /** @brief Get a free block */
  numa_block_ptr GetBlock()
  {
    numa_block_ptr block;
    if (!fBlocks.dequeue(block)) {
      block = AddBlock();
      //printf("++++++++++ NEW BLOCK %d (%d)\n", block->GetId(), block->GetNode());
    } else {
      block->Clear();
    }
    //printf("Extracted block %d (%d)\n", block->GetId(), block->GetNode());
    return block;
  }

  /** @brief Recycle a block */
  GEANT_FORCE_INLINE
  void RecycleBlock(numa_block_ptr block)
  {
//    block->Clear();
    while (!fBlocks.enqueue(block))
      ;
  }
  
  /** @brief Get a free object from the pool 
      @return A valid object reference */
  T &GetObject(size_t &index) {
    // Get the next object from the current active block
    numa_block_ptr block = CurrentBlock();
    // Hold the object
    T* obj = block->GetObject(index);
    // If the block is not yet fully distributed, return
    if (obj) return (*obj);
    // Replace distributed block
    numa_block_ptr next_free = GetBlock();
    //printf(" --- current= %d (%d) distributed, next_free= %d (%d)\n", block->GetId(), block->GetNode(), next_free->GetId(), next_free->GetNode());
    while (!fCurrent.compare_exchange_weak(block, next_free, std::memory_order_relaxed)) {
      // Retry if block is the same
      if (CurrentBlock() == block) continue;
      //printf("    no_replace for block %d (%d)\n", block->GetId(), block->GetNode());
      // Someone else replaced the block, recycle our free block
      RecycleBlock(next_free);
      break;
    }
    // Return the held object if valid, or try again with the new block
    //printf("    current is now: %d (%d)\n", CurrentBlock()->GetId(), CurrentBlock()->GetNode());
    block = CurrentBlock();
    // Blocks are large, unlikely to be emptied righ away, but you never know...
    obj = block->GetObject(index);
    if (!obj) obj = &GetObject(index);
    // There is no link anymore to the replaced block but once released by all users
    // it will go back in the block list
    return (*obj);
  }

  /** @brief Recycle an object from a block 
      @param block Block from which the object is released
      @return Block may have been recycled */
  GEANT_FORCE_INLINE bool ReleaseObject(numa_block_ptr block) {
    if (block->ReleaseObject()) {
      fNreleased++;
      RecycleBlock(block);
      //printf("ooooo  Recycling block %d (%d) Qsize = %ld\n", block->GetId(), block->GetNode(), fBlocks.size());
      return true;
    }
    return false;
  }
  
  numa_block_ptr CurrentBlock() const { return ( fCurrent.load() ); }
  
};  
} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
