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
#include <type_traits>

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

  /** @brief  GeantObjectPool copy constructor */
  GeantObjectPool(const GeantObjectPool &);

  /** @brief Assignment operator */
  GeantObjectPool &operator=(const GeantObjectPool &);

public:

  /**
   * @brief GeantObjectPool constructor
   * 
   * @param ninitial  Initial number of objects to be allocated
   * @param refobj Reference object to serve as blueprint for new allocations
   */
  GeantObjectPool(Int_t ninitial, const T *refobj = 0);

  /** @brief GeantObjectPool destructor */
  ~GeantObjectPool();
  
  /** @brief Borrow an object from the pool */
  T *Borrow();

  /**
   * @brief Create object and push it in queue
   * 
   * @param nobj Number of objects
   */
  void CreateAndPush(Int_t nobj);

  /**
   * @brief Returns back an object to the pool
   * 
   * @param obj Object that should be returned
   */
  void Return(T *obj);

  /**
   * @brief Function to set the blueprint
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
  static_assert(std::is_copy_constructible<T>::value, "Type used in GeantObjectPool must be copy constructible");
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
 * @details Returns a borrowed object to the pool.
 */
template <class T> void GeantObjectPool<T>::Return(T *obj) {
  fPool.push(obj);
}

/**
 * @details Set the "blueprint" object for the pool. If this is called, every allocation
 *  done by the pool will use the copy constructor, which is mandatory.
 */
template <class T> void GeantObjectPool<T>::SetBlueprint(const T &refobj) {
  fBlueprint = new T(refobj);
}

#endif
