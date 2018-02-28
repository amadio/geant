//===--- FactoryStore.h - GeantV --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file FactoryStore.h
 * @brief Implementation of store keeping factories one factory per user type
 * in GeantV prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_FACTORYSTORE
#define GEANT_FACTORYSTORE

#include <typeinfo>
#include <mutex>
#include <cstring>

#ifndef GEANT_FACTORY
#include "Geant/Factory.h"
#endif

/**
 * @brief Class FactoryStore
 * @details  Store of factories holding one factory per user type.
 */
class FactoryStore {
private:
  int fNclients;                   /** Number of thread clients */
  int fNFactories;                 /** Number of factories stored */
  int fCapacity;                   /** Store capacity */
  const void **fTypes;             /** Factory types */
  void **fFactories;               /** Stored factories */
  std::mutex fMutex;               /** Mutex for factory creation */
  static FactoryStore *fgInstance; /** Static instance of the store */

  /**
   * @brief Constructor for FactoryStore
   *
   * @param nclients Number of thread clients
   */
  FactoryStore(int nclients);

  /**
   * @brief Function for removing a factory
   *
   * @param islot Factory slot
   */
  void RemoveFactory(int islot);

  /**
   * @brief Constructor FactoryStore */
  FactoryStore(const FactoryStore &);

  /**
   * @brief Operator= not allowed */
  FactoryStore &operator=(const FactoryStore &);

public:
  /**
  * @brief Destructor FactoryStore */
  ~FactoryStore();

  /**
   * @brief  Function that creates one instance of the factory store
   *
   * @param nclients Number of thread clients (by default 1)
   */
  static FactoryStore *Instance(int nclients = 1);

  /**
   * @brief Templated function that return factory object
   *
   * @tparam Data type for the factory
   * @param blocksize Size of block
   */
  template <class T>
  Factory<T> *GetFactory(int blocksize, int nthreads);

  /** @brief Function that provides deletion of factory of the provided type
   *
   * @tparam Data type for the factory
  */
  template <class T>
  void DeleteFactory();
};

/**
 * @details Returns existing factory for the user type or create new one.
 * @return Factory object
 */
template <class T>
Factory<T> *FactoryStore::GetFactory(int blocksize, int nthreads)
{
  const std::type_info *type = &typeid(T);
  for (int i = 0; i < fNFactories; i++) {
    if ((const std::type_info *)fTypes[i] == type) return (Factory<T> *)fFactories[i];
  }
  Factory<T> *factory = new Factory<T>(fNclients, blocksize, nthreads);
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
  fTypes[fNFactories]       = (const void *)type;
  fFactories[fNFactories++] = factory;
  fMutex.unlock();
  return factory;
}

/**
 * @details Checks the typeid of the provided type and deletes the
 * corresponding factory if it finds one
 */
template <class T>
void FactoryStore::DeleteFactory()
{
  const std::type_info *type = &typeid(T);
  for (int i = 0; i < fNFactories; i++) {
    if ((const std::type_info *)fTypes[i] == type) {
      Factory<T> *factory = (Factory<T> *)fFactories[i];
      delete factory;
      RemoveFactory(i);
      return;
    }
  }
}
#endif
