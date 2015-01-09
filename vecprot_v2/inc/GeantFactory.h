//===--- GeantFactory.h - Geant-V -------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantFactory.h
 * @brief Implementation of factory of user objects in Geant-V prototype
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_FACTORY
#define GEANT_FACTORY

#include <Rtypes.h>
#include <vector>
#include <typeinfo>
#include <TMutex.h>
#include "dcqueue.h"
#include "WorkloadManager.h"
#include "TGeoManager.h"

using namespace std;

/**
 * @brief GeantBlock class 
 * @details  Fixed-size block of user objects stored contiguously in memory.
 * The class is used to efficiently manage user data in a multithreaded environment.
 * @tparam GeantBlock objects
*/
template <class T> class GeantBlock {
private:
  Int_t fSize;      /** Fixed size */
  Int_t fNext;      /** Number of objects in use */
  vector<T> fBlock; /** block of user objects */

public:

  /** 
   * @brief GeantBlock parameterized constructor
   * 
   * @param size  Fixed size for GeantBlock
   */
  GeantBlock(Int_t size);

  /**
   * @brief Templated destructor for GeantBlock */
  ~GeantBlock();
  
  /**
   * @brief Templated function of addition objects
   * 
   * @param p Pointer to object
   * @param index Given index to copy for
   */
  void Add(void *p, Int_t index = -1);

  /**
   * @brief ?????
   * 
   * @param index Given index for operations
   */
  T *At(Int_t index) const { return &fBlock[index]; }

  /** @brief Clear function */
  void Clear();

  /** @brief Function that check fullness of block */
  Bool_t IsFull() const { return (fNext == fSize); }

  /** @brief Function that check next free block */
  T *NextFree() { return (fNext < fSize) ? &fBlock[fNext++] : 0; }

  /** @brief Function that return size of block */
  Int_t Size() const { return fSize; }
};

/**
 * @tparam GeantBlock objects
 * @todo Some small details
 */
template <class T> GeantBlock<T>::GeantBlock(Int_t size) : fSize(size), fNext(0), fBlock() {
  fBlock.reserve(size);
  for (Int_t i = 0; i < fSize; i++)
    fBlock.push_back(T());
}

/**
 * @tparam GeantBlock objects 
 * @todo  More details?
 */
template <class T> GeantBlock<T>::~GeantBlock() {
  fBlock.clear();
}

/**
 * @details Copy the content of object pointed by p at a given index or at the next free slot.
 * @tparam GeantBlock objects  
 */
template <class T> void GeantBlock<T>::Add(void *p, Int_t index) {
  if (index < 0)
    index = fNext++;
  T &current = fBlock[index];
  current = *(T *)p;
}

/** 
 * @details Clear all data and free the block. Note that the objects are not deleted,
 * but only filled using the assignment operator from a default dummy object.
 * @tparam GeantBlock objects 
 */
template <class T> void GeantBlock<T>::Clear() {
  static T dummy;
  for (Int_t i = 0; i < fSize; i++)
    fBlock[i] = dummy;
  fNext = 0;
}

/**
 * @brief Class GeantBlockArray 
 * @details An array of blocks of user objects of the same category. Used as internal utility by factories.
 * The array access functions are critical sections.
 * 
 * @tparam Objects of GeantBlockArray type
 */
template <class T> class GeantBlockArray {
private:
  Int_t fNthreads;         /** Number of threads */
  Int_t fBlockSize;        /** Block size */
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
  GeantBlockArray(Int_t nthreads, Int_t blocksize);

  /** @brief GeantBlockArray destructor */
  ~GeantBlockArray();
  
  /**
   * @brief Operator[] 
   * 
   * @todo  Description of function (doxygen)
   * @tparam Objects of GeantBlockArray type
   * @param i ?
   */
  GeantBlock<T> *operator[](Int_t i) { return fBlocks[i]; }
  
  /**
   * @brief At Function
   * 
   * @todo  Description of function (doxygen)
   * @tparam Objects of GeantBlockArray type
   * @param i ?
   */
  GeantBlock<T> *At(Int_t i) { return fBlocks[i]; }
  
  /**
   * @brief AddAt function
   * 
   * @todo  Description of function (doxygen)
   * @param tid Track ID
   * @param block GeantBlock
   */
  void AddAt(Int_t tid, GeantBlock<T> *block) { fBlocks[tid] = block; }
};

/**
 * 
 * @tparam Objects of GeantBlockArray type
 */
template <class T>
GeantBlockArray<T>::GeantBlockArray(Int_t nthreads, Int_t blocksize)
    : fNthreads(nthreads), fBlockSize(blocksize), fBlocks(0) {
  fBlocks = new GeantBlock<T> *[nthreads];
  for (Int_t i = 0; i < nthreads; i++)
    fBlocks[i] = new GeantBlock<T>(blocksize);
}

/**
 * @tparam Objects of GeantBlockArray type
 */
template <class T> GeantBlockArray<T>::~GeantBlockArray() {
  for (Int_t i = 0; i < fNthreads; i++)
    delete fBlocks[i];
  delete[] fBlocks;
}

/**
 * @brief Class GeantFactory
 * @details Templated factory of user objects, allocated in contiguous 
 * blocks. It can serve a number of concurrent clients with id's from 0 to N.
 */
template <class T> class GeantFactory {
  friend class GeantFactoryStore;
  typedef void (*ProcessHitFunc_t)(const vector<T> &, int);

private:

  /**
   * @brief GeantFactory constructor
   * 
   * @param nthreads Number of threads
   * @param blocksize Block size
   * @param callback Callback (by default = 0)
   */
  GeantFactory(Int_t nthreads, Int_t blocksize, ProcessHitFunc_t callback = 0);
  
  /** @brief Copy constructor GeantFactory */
  GeantFactory(const GeantFactory &);
  
  /*d* @brief Operator = */
  GeantFactory &operator=(const GeantFactory &);

public:
  Int_t fNslots;                     /** Number of event slots */
  Int_t fNthreads;                   /** Max number of threads accessing the structure */
  Int_t fBlockSize;                  /** Block size */
  ProcessHitFunc_t fCallback;        /** User function to call back */
  GeantBlockArray<T> **fBlockA;      /** [fNslots] arrays of data blocks */
  dcqueue<GeantBlock<T> *> fPool;    /** pool of empty/recycled blocks */
  dcqueue<GeantBlock<T> *> fOutputs; /** Pool of filled blocks */
  
  /** GeantFactory destructor */
  ~GeantFactory();
  
  /**
   * @brief Function that add new blocks to the factory
   * 
   * @param nblocks Number of blocks
   */
  void AddFreeBlocks(Int_t nblocks);
  
  /**
   * @brief Function of freeing next block
   * 
   * @param client Client
   */
  T *NextFree(Int_t client);

  /**
   * @brief Recycle function
   * 
   * @param block Block that should be recycled
   */
  void Recycle(GeantBlock<T> *block) { fPool.push(block); }
};

/** @tparam Objects of GeantFactory type */
template <class T>
GeantFactory<T>::GeantFactory(Int_t nslots, Int_t blocksize, ProcessHitFunc_t callback)
    : fNslots(nslots), fNthreads(1), fBlockSize(blocksize), fCallback(callback), fBlockA(0),
      fPool(), fOutputs() {
  // Reserve the space for the block arrays on event slots
  fBlockA = new GeantBlockArray<T> *[fNslots];
  // Check max number of threads
  fNthreads = WorkloadManager::Instance()->GetNthreads();
  // Add 2*nclients free blocks
  AddFreeBlocks(2 * fNthreads); // why 2 ?
  for (Int_t iev = 0; iev < fNslots; iev++) {
    // One block array per slot
    fBlockA[iev] = new GeantBlockArray<T>(fNthreads, blocksize);
  }
}

/** @tparam Objects of GeantFactory type */
template <class T> GeantFactory<T>::~GeantFactory() {
  for (Int_t iev = 0; iev < fNslots; iev++)
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

/** @tparam Objects of GeantFactory type  */
template <class T> void GeantFactory<T>::AddFreeBlocks(Int_t nblocks) {
  for (Int_t i = 0; i < nblocks; i++)
    fPool.push(new GeantBlock<T>(fBlockSize));
}


/** @tparam Objects of GeantFactory type*/
template <class T> T *GeantFactory<T>::NextFree(Int_t slot) {
  // If the block is full put it
  Int_t tid = TGeoManager::ThreadId(); // maybe put in calling sequence
  GeantBlock<T> *block;
  if (fBlockA[slot]->At(tid)->IsFull()) {
    // The last entry in the block was used and filled (by the same thread)
    fOutputs.push(fBlockA[slot]->At(tid));
    fPool.wait_and_pop(block);
    fBlockA[slot]->AddAt(tid, block);
    // Keep the pool full
    if (fPool.size_async() < size_t(fNthreads))
      AddFreeBlocks(fNthreads);
  }
  return fBlockA[slot]->At(tid)->NextFree();
}
#endif
