#ifndef GEANT_FACTORY
#define GEANT_FACTORY

#include <Rtypes.h>
#include <vector>
#include <typeinfo>
#include <TMutex.h>
#include "sync_objects.h"
#include "WorkloadManager.h"
#include "TGeoManager.h"

//______________________________________________________________________________
// GeantBlock - a fixed-size block of user objects stored contiguously in memory
//                 The class is used to efficiently manage user data in a
//                 multithreaded environment.
//______________________________________________________________________________

using namespace std;

template <class T>
class GeantBlock {
private:
   Int_t                 fSize;      // fixed size
   Int_t                 fNext;      // number of objects in use
   vector<T>             fBlock;     // block of user objects

public:
   GeantBlock(Int_t size);
   ~GeantBlock();

   void          Add(void *p, Int_t index=-1);
   T            *At(Int_t index) const {return &fBlock[index];}
   void          Clear();
   Bool_t        IsFull() const {return (fNext==fSize);}
   T            *NextFree() {return (fNext<fSize) ? &fBlock[fNext++] : 0;}
   Int_t         Size() const {return fSize;}
};

template <class T>
GeantBlock<T>::GeantBlock(Int_t size)
              :fSize(size),
               fNext(0),
               fBlock()
{
// Constructor.
   fBlock.reserve(size);
   for (Int_t i=0; i<fSize; i++)
      fBlock.push_back(T());
}

template <class T>
GeantBlock<T>::~GeantBlock()
{
// Destructor
   fBlock.clear();
}   

template <class T>
void GeantBlock<T>::Add(void *p, Int_t index)
{
// Copy the content of object pointed by p at a given index or at the next
// free slot.
   if (index<0) index = fNext++;
   T &current = fBlock[index];
   current = *(T*)p;
}
   
template <class T>
void GeantBlock<T>::Clear()
{
// Clear all data and free the block. Note that the objects are not deleted,
// but only filled using the assignment operator from a default dummy object.
   static T dummy;
   for (Int_t i=0; i<fSize; i++) fBlock[i] = dummy;
   fNext = 0;
}

//______________________________________________________________________________
// GeantBlockArray - An array of blocks of user objects of the same category.
//                   Used as internal utility by factories.
//                   The array access functions are critical sections.
//______________________________________________________________________________

template <class T>
class GeantBlockArray {
private:
   Int_t                  fNthreads;  // Number of threads
   Int_t                  fBlockSize; // Block size
   GeantBlock<T>        **fBlocks;    // array of blocks used by different threads

   GeantBlockArray(const GeantBlockArray&);
   GeantBlockArray& operator=(const GeantBlockArray&);
public:
   GeantBlockArray(Int_t nthreads, Int_t blocksize);
   ~GeantBlockArray();
   
   GeantBlock<T> *operator[](Int_t i) {return fBlocks[i];}
   GeantBlock<T> *At(Int_t i) {return fBlocks[i];}
   void           AddAt(Int_t tid, GeantBlock<T> *block) {fBlocks[tid] = block;}
};   

template <class T>
GeantBlockArray<T>::GeantBlockArray(Int_t nthreads, Int_t blocksize)
                   :fNthreads(nthreads), fBlockSize(blocksize), fBlocks(0)
{
// Constructor
   fBlocks = new GeantBlock<T>*[nthreads];
   for (Int_t i=0; i<nthreads; i++) fBlocks[i] = new GeantBlock<T>(blocksize);
}   
   
template <class T>
GeantBlockArray<T>::~GeantBlockArray()
{
// Destructor
   for (Int_t i=0; i<fNthreads; i++) delete fBlocks[i];
   delete [] fBlocks;
}   

   
//______________________________________________________________________________
// GeantFactory - Templated factory of user objects, allocated in contiguous 
//                blocks. It can serve a number of concurrent clients with
//                id's from 0 to N.
//______________________________________________________________________________

template <class T>
class GeantFactory {
   friend class GeantFactoryStore;
   typedef void (*ProcessHitFunc_t)(const vector<T>&, int);
private:
   GeantFactory(Int_t nthreads, Int_t blocksize, ProcessHitFunc_t callback=0);
   GeantFactory(const GeantFactory &);
   GeantFactory &operator=(const GeantFactory &);
public:
   Int_t                  fNslots;      // Number of event slots
   Int_t                  fNthreads;    // Max number of threads accessing the structure
   Int_t                  fBlockSize;   // Block size
   ProcessHitFunc_t       fCallback;    // User function to call back
   GeantBlockArray<T>   **fBlockA;      //[fNslots] arrays of data blocks
   dcqueue< GeantBlock<T> > fPool;        // pool of empty/recycled blocks
   dcqueue< GeantBlock<T> > fOutputs;     // pool of filled blocks

    ~GeantFactory();
   
   void          AddFreeBlocks(Int_t nblocks);
   T            *NextFree(Int_t client);
   void          Recycle(GeantBlock<T>* block) {fPool.push(block);}
};

template <class T>
GeantFactory<T>::GeantFactory(Int_t nslots, Int_t blocksize, ProcessHitFunc_t callback)
              :fNslots(nslots),
               fNthreads(1),
               fBlockSize(blocksize),
               fCallback(callback),
               fBlockA(0),
               fPool(),
               fOutputs()
{
// Constructor.
   // Reserve the space for the block arrays on event slots
   fBlockA = new GeantBlockArray<T>*[fNslots];
   // Check max number of threads
   fNthreads = WorkloadManager::Instance()->GetNthreads();
   // Add 2*nclients free blocks
   AddFreeBlocks(2*fNthreads); // why 2 ?
   for (Int_t iev=0; iev<fNslots; iev++) {
      // one block array per slot
      fBlockA[iev] = new GeantBlockArray<T>(fNthreads, blocksize);
   }
}

template <class T>
GeantFactory<T>::~GeantFactory()
{
// Destructor
   for (Int_t iev=0; iev<fNslots; iev++) delete [] fBlockA[iev];
   delete [] fBlockA;
   while (!fPool.empty())   {delete fPool.back(); fPool.pop_back();}
   while (!fOutputs.empty())   {delete fOutputs.back(); fOutputs.pop_back();}
}

template <class T>
void GeantFactory<T>::AddFreeBlocks(Int_t nblocks)
{
// Add new blocks to the factory
   for (Int_t i=0; i<nblocks; i++) fPool.push(new GeantBlock<T>(fBlockSize));
}   

template <class T>
T *GeantFactory<T>::NextFree(Int_t slot)
{
// Returns the next free object
   // If the block is full put it 
   Int_t tid = TGeoManager::ThreadId(); // maybe put in calling sequence
   if (fBlockA[slot]->At(tid)->IsFull()) {
      // The last entry in the block was used AND filled (by the same thread)
      fOutputs.push(fBlockA[slot]->At(tid));
      fBlockA[slot]->AddAt(tid, fPool.wait_and_pop());
      // Keep the pool full
      if (fPool.size_async() < fNthreads) AddFreeBlocks(fNthreads);
   }   
   return fBlockA[slot]->At(tid)->NextFree();
}
#endif
