//===--- GeantObjectPool.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantObjectPool.h
 * @brief Concurrent pool of generic pre-alocated objects providing the borrow/return functionality
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_OBJECTPOOL
#define GEANT_OBJECTPOOL
#include "priority_queue.h"

/**
 * @brief GeantObjectPool class
 * @details Concurrent pool of generic pre-alocated objects providing the borrow/return 
 * functionality
 */
template <class T> class GeantObjectPool {
  Geant::priority_queue<T> fPool; /** Concurrent queue used to pool objects */
  T *fBlueprint;                  /** Blueprint object from which each new allocation */
                                  /** will be copied from. Requires working CC. */
private:

  /** @brief  GeantObjectPool constructor */
  GeantObjectPool(const GeantObjectPool &);

  /** @brief operator= */
  GeantObjectPool &operator=(const GeantObjectPool &);

public:

  /**
   * @brief GeantObjectPool constructor
   * 
   * @param ninitial  Initial number ?
   * @param refobj Referenced object
   */
  GeantObjectPool(Int_t ninitial, const T *refobj = 0);

  /** @brief GeantObjectPool destructor */
  ~GeantObjectPool();
  
  /** @brief Function of borrowing */
  T *Borrow();

  /**
   * @brief Create object and push it in queue
   * 
   * @param nobj Number of objects
   */
  void CreateAndPush(Int_t nobj);

  /**
   * @brief Return function
   * 
   * @param obj Object that should be returned
   */
  void Return(T *obj);

  /**
   * @brief Function that set blueprint
   * 
   * @param refobj Referenced object
   */
  void SetBlueprint(const T &refobj);
};

/**
 * @details It allows to define the initial capacity and to provide
 * a blueprint for the pooled objects.
 */
template <class T>
GeantObjectPool<T>::GeantObjectPool(Int_t ninitial, const T *refobj)
    : fPool(), fBlueprint(0) {
  //   if (!refobj) fBlueprint = new T();  // assumes default ctor
  //   else
  fBlueprint = new T(*refobj); // uses CC
  CreateAndPush(ninitial);
}

/** 
 * @details GeantObjectPool destructor.
 * Calls also the destructor of remaining objects 
 */
template <class T> GeantObjectPool<T>::~GeantObjectPool() {
  fPool.delete_content();
}

/**
 * @brief Create object and push it in queue
 * @details Create nobjects and push them in the queue. This should be done only at
 * initialization
 */
template <class T> void GeantObjectPool<T>::CreateAndPush(Int_t nobj) {
  for (Int_t i = 0; i < nobj; i++)
    fPool.push(new T(*fBlueprint));
}

/**
 * @details Borrow an object from the pool.
 * @return Borrowed object
 */
template <class T> T *GeantObjectPool<T>::Borrow() {
  T *obj = fPool.try_pop();
  if (!obj)
    obj = new T(*fBlueprint);
  return obj;
}

/**
 * @details Returns a borrowed object.
 */
template <class T> void GeantObjectPool<T>::Return(T *obj) {
  fPool.push(obj);
}

/**
 * @details Set the "blueprint" object for the pool. If this is called, every allocation
 *  done by the pool will use the copy constructor, else it will fallback on
 *  calling the default constructor, which is in this case mandatory.
 */
template <class T> void GeantObjectPool<T>::SetBlueprint(const T &refobj) {
  fBlueprint = new T(refobj);
}

#endif
