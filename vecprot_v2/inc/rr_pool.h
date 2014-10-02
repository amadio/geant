#ifndef GEANT_RR_POOL
#define GEANT_RR_POOL
#include <vector>
#include <iostream>
#if __cplusplus >= 201103L
#include <atomic>
#endif

//!_____________________________________________________________________________
//! Generic atomic round robin pool
//!_____________________________________________________________________________

/**
   Pool of arbitrary objects allocated on the heap. The pool is created by 
   providing a number of slots which will be selected in a round robin manner.
   The user must create a blueprint object before creating the pool and provide 
   it to the pool constructor. Whenever the selected slot has no stored object, 
   one will be created using the copy constructor of the generic type. The 
   capacity represents the number of object for each slot which will be created 
   at construction time.
   
   For an optimal usage in a multithreaded environment, one should request a 
   number of slots greated than the number of user threads. The capacity should 
   be greater than the number of objects from the pool held at any moment by 
   any thread; this will minimize run time allocations. Any object from the pool 
   has to be acquired using the borrow() method and released using the return() 
   method.   
*/

template <class T>
class rr_pool {
public:
   typedef std::vector<T*> rrlist_t;
   int                fNslots;       //!< Number of slots
   int                fCapacity;     //!< Nomber of object to be preallocated per slot
   T                 *fBlueprint;    //!< Blueprint object to be copied for allocations 
#if __cplusplus >= 201103L
   std::atomic_int    fRRcount;      //!< Round robin counter
   std::atomic_int    fNcalls;       //!< Number of calls to borrow/release
   std::atomic_int    fInew;         //!< Number of allocations done (overhead)
   std::atomic_int    fLockHit;      //!< Number of lock hits (overhead)
   std::vector<std::atomic_flag*>   fLocks;  //!< Memory barriers per slot
#endif
   std::vector<rrlist_t*>  fLists;  //!< List of lists of objects per slot

private:
   rr_pool(const rr_pool&);
   rr_pool &operator=(const rr_pool&);

public:
   rr_pool(int nslots, int capacity, const T *refobj);
   ~rr_pool();
   
   T*                 borrow();
   void               release(T* obj);
   int                next_slot();
   void               statistics();
};

// *** Implementation ***
//______________________________________________________________________________
template <class T>
rr_pool<T>::rr_pool(int nslots, int capacity, const T *refobj) 
           :fNslots(nslots), fCapacity(capacity), fBlueprint(0), fRRcount(0), 
            fNcalls(0), fInew(0), fLockHit(0), fLocks(), fLists()
{
//! Constructor.
   rrlist_t *list;
   fBlueprint = new T(*refobj); // uses CC
   for (auto i=0; i<nslots; ++i) {
      list = new rrlist_t;
      list->reserve(capacity);
      for (auto j=0; j<capacity; ++j) list->push_back(new T(*fBlueprint));
      fLists.push_back(list);
      std::atomic_flag *lock = new std::atomic_flag;
      // Must expliticit initialize the atomic_flag to avoid valgring errors.
      lock->clear();
      fLocks.push_back(lock);
   }
}

//______________________________________________________________________________
template <class T>
rr_pool<T>::~rr_pool()
{
//! Destructor.
   rrlist_t *list;
   for (auto it = fLists.begin(); it !=  fLists.end(); ++it) {
      list = *it;
      for (auto it1 = list->begin(); it1 != list->end(); ++it1) delete (*it1);
      list->clear();
      delete list;
   }
   for (auto it = fLocks.begin(); it != fLocks.end(); ++it) delete (*it);
   fLocks.clear();
   fLists.clear();
}

//______________________________________________________________________________
template <class T>
T* rr_pool<T>::borrow()
{
//! Borrows an object from the pool. If none available create one.
   int islot = next_slot();
   fNcalls++;
   while (fLocks[islot]->test_and_set(std::memory_order_acquire)) {
      fLockHit++;
      islot = next_slot();
   }   
   rrlist_t *list = fLists[islot];
   T *retobj = 0;
   if (list->empty()) {
      fLocks[islot]->clear(std::memory_order_release);
      fInew++;
      return new T(*fBlueprint);
   }
   retobj = list->back();
   list->pop_back();
   fLocks[islot]->clear(std::memory_order_release);
   return retobj;
}

//______________________________________________________________________________
template <class T>
void rr_pool<T>::release(T *obj)
{
//! Returns an object to the pool.
   int islot = next_slot();
   fNcalls++;
   while (fLocks[islot]->test_and_set(std::memory_order_acquire)) {
      fLockHit++;
      islot = next_slot();
   }   
   rrlist_t *list = fLists[islot];
   list->push_back(obj);   
   fLocks[islot]->clear(std::memory_order_release);
}

//______________________________________________________________________________
template <class T>
int rr_pool<T>::next_slot()
{
//! Returns next slot in round robin.
   int nslots = fNslots;
   int slot = fRRcount++;
   fRRcount.compare_exchange_strong(nslots, 0);
   return (slot % fNslots);
}

//______________________________________________________________________________
template <class T>
void rr_pool<T>::statistics()
{
//! Print usage statistics
   std::cout << "=== # RR slots        : " << fNslots << std::endl;
   std::cout << "=== # pre-allocs/slot : " << fCapacity << std::endl;
   std::cout << "=== # pool calls      : " << fNcalls.load() << std::endl;
   std::cout << "=== # allocations     : " << fInew.load() << std::endl;
   std::cout << "=== # locked spins    : " << fLockHit.load() << std::endl;
}   

#endif
