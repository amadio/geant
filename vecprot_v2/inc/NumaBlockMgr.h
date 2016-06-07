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

template <typename T> class NumaBlockMgr {
  using size_t = std::size_t size_t;
  using numa_block_ptr = NumaBlock<T>*;
  using queue_t = mpmc_bounded_queue<numa_block_ptr>;

  static size_t const cacheline_size = 64;
  static size_t const queue_buff_size = 4096;
  typedef char cacheline_pad_t[cacheline_size];

private:
  std::atomic_flag fLock; // Atomic lock

  cacheline_pad_t pad0_;   //! Padding to protect the other data from the hot cache line above

  std::atomic<numa_block_ptr> fCurrent; // Current block being distributed

  cacheline_pad_t pad1_;   //! Padding to protect the other data from the hot cache line above
  
  int             fNode;  // Numa node id
  size_t          fBlockSize;  // Numa block size
  queue_t fBlocks;  // Queue of free blocks

public:

  /** @brief Constructor providing number of blocks to be initially created */
  NumaBlockMgr(int numa_node, size_t bsize, int ninitial) 
    : fLock(), fCurrent(nullptr), fNode(numa_node), fBlockSize(bsize), 
      fBlocks(queue_buff_size), fUsed(queue_buff_size)
  {
    fLock.clear();
    assert(ninitial <= queue_buff_size);
    for (auto i=0; i<ninitial; ++i) {
      fBlocks.enqueue(AddBlock());
    }
  }
  
  /** @brief Add a free block */
  numa_block_ptr AddBlock()
  {
    return ( NumaBlock<T>::MakeInstance(fBlockSize, fNode) );
  }

  /** @brief Get a free block */
  numa_block_ptr GetBlock()
  {
    numa_block_ptr block;
    if (!fBlocks.dequeue(block))
      block = AddBlock();
    return block;
  }

  
  /** @brief Get a free object from the pool 
      @return A valid object reference */
  T & const GetObject() {
    // Get the next object from the current active block
    NumaBlock_ptr block = CurrentBlock();
    T* obj = block->GetObject();
    if (obj) return (*obj);
    // Block fully distributed, replace with free one
    NumaBlock_ptr next_free = GetBlock();
    while (!fCurrent.compare_exchange_weak(block, next_free, std::memory_order_relaxed)) {
      // If not managing to replace the old block, someone else may have done it
      block = CurrentBlock();
      obj = block->GetObject();
      if (obj) {
        // Someone else indeed replaced the block, recycle our free block
        assert(fBlocks.enqueue(next_free));
        return (*obj);
      }
    }
    // Blocks are large, unlikely to be emptied righ away, but you never know...
    assert( (obj = block->GetObject()) );
    // There is no link anymore to the replaced block but once released by all users
    // it will go back in the empty list
    assert(fUsed.enqueue(block));
    return obj;
  }

  /** @brief Recycle an object from a block 
      @param block Block from which the object is released
      @return Block may have been recycled */
  GEANT_INLINE bool ReleaseObject(NumaBlock_ptr block) { 
    if (block->ReleaseObject() == 0) {
      fBlocks.enqueue(block);
      return true;
    }
    return false;
  }
  
  GEANT_INLINE void Lock() { while (fLock.test_and_set(std::memory_order_acquire)) ; }
  GEANT_INLINE void Unlock() { fLock.clear(std::memory_order_release); }
  NumaBlock_ptr CurrentBlock() const { return ( fCurrent.load() ); }
  
};  
} // Geant
