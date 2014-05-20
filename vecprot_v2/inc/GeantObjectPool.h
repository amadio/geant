#ifndef GEANT_OBJECTPOOL
#define GEANT_OBJECTPOOL
#include "sync_objects.h"

//______________________________________________________________________________
// Concurrent pool of generic pre-alocated objects providing the borrow/return
// functionality
//______________________________________________________________________________

template <class T>
class GeantObjectPool {
   dcqueue<T>        fPool;      // Concurrent queue used to pool objects
public:
   GeantObjectPool(Int_t ninitial) : fPool() {CreateAndPush(ninitial);}
   ~GeantObjectPool();
   
   T          *Borrow();
   void        CreateAndPush(Int_t nobj);
   void        Return(T *obj);
};
   
//______________________________________________________________________________
template <class T>
GeantObjectPool<T>::~GeantObjectPool()
{
// destructor. Calls also the destructor of remaining objects
   fPool.delete_content();
}

//______________________________________________________________________________
template <class T>
void GeantObjectPool<T>::CreateAndPush(Int_t nobj)
{
// Create nobjects and push them in the queue. This should be done only at
// initialization
   for (Int_t i=0; i<nobj; i++) fPool.push(new T());
}

//______________________________________________________________________________
template <class T>
 T* GeantObjectPool<T>::Borrow()
{
// Borrow an object from the pool.
   T *obj = fPool.try_pop();
   if (!obj) obj = new T();
   return obj;
}    

//______________________________________________________________________________
template <class T>
 void GeantObjectPool<T>::Return(T *obj)
{
// Return a borrowed object.
   fPool->push(obj);
}   
#endif
