//===--- GeantFactoryStore.h - GeantV --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantFactoryStore.h
 * @brief Implementation of store keeping factories one factory per user type 
 * in GeantV prototype
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_FACTORYSTORE
#define GEANT_FACTORYSTORE

#include <typeinfo>
#include <mutex>

#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

/** 
 * @brief Class GeantFactoryStore
 * @details  Store of factories holding one factory per user type.
 */
class GeantFactoryStore {
private:
  int          fNclients;               /** Number of thread clients */
  int          fNFactories;             /** Number of factories stored */
  int          fCapacity;               /** Store capacity */
  const void **fTypes;                  /** Factory types */
  void       **fFactories;              /** Stored factories */
  std::mutex   fMutex;                  /** Mutex for factory creation */
  static GeantFactoryStore *fgInstance; /** Static instance of the store */
  
  /**
   * @brief Constructor for GeantFactoryStore
   * 
   * @param nclients Number of thread clients 
   */
  GeantFactoryStore(int nclients);
  
  /**
   * @brief Function for removing a factory
   * 
   * @param islot Factory slot
   */
  void RemoveFactory(int islot);

  /**
   * @brief Constructor GeantFactoryStore */
  GeantFactoryStore(const GeantFactoryStore &);

  /**
   * @brief Operator= not allowed */
  GeantFactoryStore &operator=(const GeantFactoryStore &);
public:
  
  /** 
  * @brief Destructor GeantFactoryStore */
  ~GeantFactoryStore();
  
  /**
   * @brief  Function that creates one instance of the factory store
   * 
   * @param nclients Number of thread clients (by default 1)
   */
  static GeantFactoryStore *Instance(int nclients = 1);
  
  /**
   * @brief Templated function that return factory object
   * 
   * @tparam Data type for the factory
   * @param blocksize Size of block
   */
  template <class T> GeantFactory<T> *GetFactory(int blocksize,WorkloadManager *wMgr);
  
  /** @brief Function that provides deletion of factory of the provided type 
   *
   * @tparam Data type for the factory
  */
  template <class T> void DeleteFactory();
};

/**
 * @details Returns existing factory for the user type or create new one.
 * @return Factory object
 */
template <class T> GeantFactory<T> *GeantFactoryStore::GetFactory(int blocksize, WorkloadManager *wMgr) {
  const std::type_info *type = &typeid(T);
  for (int i = 0; i < fNFactories; i++) {
    if ((const std::type_info *)fTypes[i] == type)
      return (GeantFactory<T> *)fFactories[i];
  }
  GeantFactory<T> *factory = new GeantFactory<T>(fNclients, blocksize, wMgr);
  fMutex.lock();
  if (fNFactories == fCapacity) {
    // Resize arrays
    const void **types = new const void *[2 * fCapacity];
    memset(types, 0, 2 * fCapacity * sizeof(void *));
    memcpy(types, fTypes, fNFactories * sizeof(void *));
    void **factories = new void *[2 * fCapacity];
    memset(factories, 0, 2 * fCapacity * sizeof(void *));
    memcpy(factories, fFactories, fNFactories * sizeof(void *));
    delete[] fTypes;
    fTypes = types;
    delete[] fFactories;
    fFactories = factories;
    fCapacity *= 2;
  }
  fTypes[fNFactories] = (const void *)type;
  fFactories[fNFactories++] = factory;
  fMutex.unlock();
  return factory;
}

/**
 * @details Checks the typeid of the provided type and deletes the
 * corresponding factory if it finds one
 */
template <class T> void GeantFactoryStore::DeleteFactory() {
  const std::type_info *type = &typeid(T);
  for (int i = 0; i < fNFactories; i++) {
    if ((const std::type_info *)fTypes[i] == type) {
      GeantFactory<T> *factory = (GeantFactory<T> *)fFactories[i];
      delete factory;
      RemoveFactory(i);
      return;
    }
  }
}
#endif
