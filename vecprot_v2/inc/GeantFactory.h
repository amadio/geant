//===--- GeantFactory.h - GeantV -------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantFactory.h
 * @brief Implementation of factory of user objects in Geant-V prototype
 * @details The file contains the template definitions of:
 * GeantBlock - a fixed-size vector of objects having user-defined type
 * GeantBlockArray - an array of GeantBlock objects (one per thread)
 * GeantFactory - a factory created on demand and handling GeantBlockArray
 *  objects for a fixed number of event slots
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_FACTORY
#define GEANT_FACTORY
#ifndef GEANTV_MIC
#include <Rtypes.h>
#endif

#include <vector>
#include <cassert>
#include <typeinfo>
#include <type_traits>
#include "dcqueue.h"
#include "WorkloadManager.h"
#ifdef USE_VECGEOM_NAVIGATOR
#else
#include "TGeoManager.h"
#endif

using namespace std;

/**
 * @brief GeantBlock class
 * @details  Fixed-size array of data blocks of user-defined type stored contiguously
 * in memory. The user type should be POD-like (no virtual table) and must have a
 * public default and copy contructors.The class is used to efficiently manage user
 * data in a multithreaded environment.
 *
 * @tparam User-defined type
*/
template <typename T> class GeantBlock {
private:
  int fSize;        /** Fixed size */
  int fNext;        /** Index of next free block */
  vector<T> fBlock; /** vector of user block objects */

public:
  /**
   * @brief GeantBlock constructor
   * @details The constructor makes sure that the provided type does not have
   * a virtual table, that it has a default constructor and that it is copiable.
   *
   * @param size  Size for the array of blocks
   */
  GeantBlock() : fSize(0), fNext(0), fBlock() {
    static_assert(!std::is_polymorphic<T>::value, "Cannot use polymorphic types as GeantBlock");
    static_assert(std::is_default_constructible<T>::value, "Type used in GeantBlock must have default ctor.");
    static_assert(std::is_copy_constructible<T>::value, "Type used in GeantBlock must be copy constructible");
  }

  /**
   * @brief Destructor for GeantBlock */
  ~GeantBlock() { fBlock.clear(); }

  /**
    * @brief Initialize the block for the given size.
    * @details Reserves the vector of user block objects and fills it with default objects of type T.
    *
    * @param size  Size for the array of blocks
    */
  void Initialize(int size) {
    assert(size > 0);
    fSize = size;
    fBlock.reserve(size);
    for (int i = 0; i < size; i++)
      fBlock.push_back(T());
  }

  /**
   * @brief Add an object of type T at a given index in the block vector.
   * @details Copy the content of object pointed by p at a given index or at the
   * next free slot if the input index is negative.
   *
   * @tparam GeantBlock object type
   *
   * @param p Pointer to object to be copied
   * @param index Index in the block vector where to copy to
   */
  void Add(T *p, int index = -1) {
    if (index < 0)
      index = fNext++;
    T &current = fBlock[index];
    current = *p;
  }

  /**
   * @brief Return the pointer to the object stored at a given block index.
   *
   * @param index Index to read from
   */
  const T *At(int index) const { return &fBlock[index]; }

  /** @brief Clear function
  * @details Clear all data and free the block. Note that the objects are not
  * deleted, but only filled using the assignment operator from a default
  * dummy object.
  *
  * @tparam GeantBlock object type
  */
  void Clear() {
    static T dummy;
    for (int i = 0; i < fSize; i++)
      fBlock[i] = dummy;
    fNext = 0;
  }

  /** @brief Function checking if the block is full */
  bool IsFull() const { return (fNext == fSize); }

  /** @brief Function getting next free block */
  T *NextFree() { return (fNext < fSize) ? &fBlock[fNext++] : 0; }

  /** @brief Function that returns the size of block */
  int Size() const { return fSize; }
};

/**
 * @brief Class GeantBlockArray
 * @details An array of blocks of user objects of the same category. Used as internal utility by factories.
 * The array access functions are critical sections.
 *
 * @tparam Object type to be stored in blocks
 */
template <typename T> class GeantBlockArray {
private:
  int fNthreads;           /** Number of threads */
  int fBlockSize;          /** Block size */
  GeantBlock<T> **fBlocks; /** Array of blocks used by different threads */

  /**
   * @brief Copy constructor GeantBlockArray
   *
   * @todo Still needs to be implemented
   */
  GeantBlockArray(const GeantBlockArray &);

  /**
  * @brief Implementation of operator=
  *
  * @todo Still needs to be implemented
  */
  GeantBlockArray &operator=(const GeantBlockArray &);

public:
  /**
   * @brief GeantBlockArray parametrized constructor
   *
   * @param nthreads Number of threads
   * @param blocksize Block size
   */
  GeantBlockArray(int nthreads, int blocksize) : fNthreads(nthreads), fBlockSize(blocksize), fBlocks(0) {
    fBlocks = new GeantBlock<T> *[nthreads];
    for (int i = 0; i < nthreads; i++) {
      fBlocks[i] = new GeantBlock<T>();
      fBlocks[i]->Initialize(blocksize);
    }
  }

  /** @brief GeantBlockArray destructor */
  ~GeantBlockArray() {
    for (int i = 0; i < fNthreads; i++)
      delete fBlocks[i];
    delete[] fBlocks;
  }

  /**
   * @brief Operator[]
   *
   * @tparam Object type to be stored in blocks
   * @param i Index to be accessed
   * @return Pointer to block
   */
  GeantBlock<T> *operator[](int i) { return fBlocks[i]; }

  /**
   * @brief Read block at a given index
   *
   * @tparam Object type to be stored in blocks
   * @param i Index to be accessed
   * @return Pointer to block
   */
  GeantBlock<T> *At(int i) { return fBlocks[i]; }

  /**
   * @brief Add a block at a given index
   *
   * @tparam Object type to be stored in blocks
   * @param tid Thread id
   * @param block GeantBlock pointer
   */
  void AddAt(int tid, GeantBlock<T> *block) { fBlocks[tid] = block; }
};

/**
 * @brief Class GeantFactory
 * @details Templated factory of user objects, allocated in contiguous
 * blocks. It can serve a number of concurrent clients with id's from 0 to N.
 */
template <typename T> class GeantFactory {
  friend class GeantFactoryStore;
  typedef void (*ProcessHitFunc_t)(const vector<T> &, int);

private:
  /**
   * @brief GeantFactory parameterised constructor. Can only be called by a
   * GeantFactoryStore instance.
   *
   * @param nthreads Number of threads
   * @param blocksize Block size
   * @param callback Callback (by default = 0)
   */
  GeantFactory(int nslots, int blocksize, ProcessHitFunc_t callback = 0)
    : fNslots(nslots), fNthreads(1), fBlockSize(blocksize), fCallback(callback), fBlockA(0), fPool(), fOutputs(),
      queue_per_thread(false) {
    // Reserve the space for the block arrays on event slots
    fBlockA = new GeantBlockArray<T> *[fNslots];
    // Check max number of threads
    fNthreads = WorkloadManager::Instance()->GetNthreads();

    //
    fPoolArray = new deque<GeantBlock<T> *> [fNthreads];
    fOutputsArray = new deque<GeantBlock<T> *> [fNthreads];
    
    // Add 2*nclients free blocks (2?)
    AddFreeBlocks(2 * fNthreads);
    
    // Add 2(?) free blocks per thread
    for(int i=0; i<fNthreads; i++) AddFreeBlocks(2, i);

    for (int iev = 0; iev < fNslots; iev++) {
      // One block array per slot
      fBlockA[iev] = new GeantBlockArray<T>(fNthreads, blocksize);
    }
  }

  /** @brief Copy constructor GeantFactory */
  GeantFactory(const GeantFactory &);

  /** @brief Operator= */
  GeantFactory &operator=(const GeantFactory &);

public:
  int fNslots;                       /** Number of event slots */
  int fNthreads;                     /** Max number of threads accessing the factory */
  int fBlockSize;                    /** Block size */
  ProcessHitFunc_t fCallback;        /** User function to call back */
  GeantBlockArray<T> **fBlockA;      /** [fNslots] arrays of data blocks */
  dcqueue<GeantBlock<T> *> fPool;    /** pool of empty/recycled blocks */
  dcqueue<GeantBlock<T> *> fOutputs; /** Pool of filled blocks */
  
  bool queue_per_thread;
  std::deque<GeantBlock<T> *> *fPoolArray;    /** [fNthreads] array of queues (per thread) of empty/recycled blocks */
  std::deque<GeantBlock<T> *> *fOutputsArray; /** [fNthreads] array of queues (per thread) of filled blocks */

  /** @brief GeantFactory destructor */
  ~GeantFactory() {
    for (int iev = 0; iev < fNslots; ++iev)
      delete[] fBlockA[iev];
    delete[] fBlockA;
    while (!fPool.empty()) {
      delete fPool.back();
      fPool.pop_back();
    }
    while (!fOutputs.empty()) {
      delete fOutputs.back();
      fOutputs.pop_back();
    }
  }

  /**
   * @brief Function to add new blocks to the factory
   *
   * @param nblocks Number of blocks to be added
   */
  void AddFreeBlocks(int nblocks) {
    for (int i = 0; i < nblocks; i++) {
      GeantBlock<T> *block = new GeantBlock<T>();
      block->Initialize(fBlockSize);
      fPool.push(block);
    }
  }

  /**
   * @brief Function to add new blocks to the pool for the given thread
   *
   * @param nblocks Number of blocks to be added
   */
  void AddFreeBlocks(int nblocks, int tid) {
    for (int i = 0; i < nblocks; i++) {
      GeantBlock<T> *block = new GeantBlock<T>();
      block->Initialize(fBlockSize);
      fPoolArray[tid].push_front(block);
    }
  }
  
  
  /**
   * @brief Function for getting the next free block
   *
   * @param slot Event slot id
   * @param tid Thread id
   */
  T *NextFree(int slot, int tid) {
    GeantBlock<T> *block;
    if (fBlockA[slot]->At(tid)->IsFull()) {
      // The last entry in the block was used and filled (by the same thread)     
      if(queue_per_thread)
	{
	  fOutputsArray[tid].push_front(fBlockA[slot]->At(tid));
	  block = fPoolArray[tid].back();
	  fPoolArray[tid].pop_back();
	  fBlockA[slot]->AddAt(tid, block);
	  // Keep the pool full
 	  if (fPoolArray[tid].size() < 2) // (2?)
	    AddFreeBlocks(2, tid);
	}
      else
	{
	  fOutputs.push(fBlockA[slot]->At(tid));
	  fPool.wait_and_pop(block);
	  fBlockA[slot]->AddAt(tid, block);
	  // Keep the pool full
	  if (fPool.size_async() < size_t(fNthreads))
	    AddFreeBlocks(fNthreads);
	}
    }
    return fBlockA[slot]->At(tid)->NextFree();
  }
    
  /**
   * @brief Recycle function
   *
   * @param block Block that should be recycled
   */
  void Recycle(GeantBlock<T> *block) {
    block->Clear();
    fPool.push(block);
  }
  
  /**
   * @brief Recycle function for thread-local pool
   *
   * @param block Block that should be recycled
   */
  void Recycle(GeantBlock<T> *block, int tid) {
    block->Clear();
    fPoolArray[tid].push_front(block);
  }
 
};

#endif
