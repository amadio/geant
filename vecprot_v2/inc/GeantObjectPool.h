#ifndef GEANT_OBJECTPOOL
#define GEANT_OBJECTPOOL
#include "priority_queue.h"

//______________________________________________________________________________
// Concurrent pool of generic pre-alocated objects providing the borrow/return
// functionality
//______________________________________________________________________________

template <class T>
class GeantObjectPool {
   Geant::priority_queue<T>        fPool;      // Concurrent queue used to pool objects
   T                *fBlueprint; // Blueprint object from which each new allocation 
                                 // will be copied from. Requires working CC.
private:
   GeantObjectPool(const GeantObjectPool&);
   GeantObjectPool &operator=(const GeantObjectPool&);
public:
   GeantObjectPool(Int_t ninitial, const T *refobj=0);
   ~GeantObjectPool();
   
   T          *Borrow();
   void        CreateAndPush(Int_t nobj);
   void        Return(T *obj);
   void        SetBlueprint(const T &refobj);
};

//______________________________________________________________________________
template <class T>
GeantObjectPool<T>::GeantObjectPool(Int_t ninitial, const T *refobj)
                   :fPool(), fBlueprint(0)
{
// Constructor. It allows to define the initial capacity and to provide
// a blueprint for the pooled objects.
//   if (!refobj) fBlueprint = new T();  // assumes default ctor
//   else         
   fBlueprint = new T(*refobj); // uses CC
   CreateAndPush(ninitial);
}   
   
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
   for (Int_t i=0; i<nobj; i++) fPool.push(new T(*fBlueprint));
}

//______________________________________________________________________________
template <class T>
 T* GeantObjectPool<T>::Borrow()
{
// Borrow an object from the pool.
   T *obj = fPool.try_pop();
   if (!obj) obj = new T(*fBlueprint);
   return obj;
}    

//______________________________________________________________________________
template <class T>
 void GeantObjectPool<T>::Return(T *obj)
{
// Return a borrowed object.
   fPool.push(obj);
}   

//______________________________________________________________________________
template <class T>
void GeantObjectPool<T>::SetBlueprint(const T &refobj)
{
// Set the "blueprint" object for the pool. If this is called, every allocation
// done by the pool will use the copy constructor, else it will fallback on
// calling the default constructor, which is in this case mandatory.
   fBlueprint = new T(refobj);
}

#endif
