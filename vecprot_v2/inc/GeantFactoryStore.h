//===--- GeantFactoryStore.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantFactoryStore.h
 * @brief Implementation of store of factories holding by one factory per user type in Geant-V prototype
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_FACTORYSTORE
#define GEANT_FACTORYSTORE

#include <Rtypes.h>
#include <typeinfo>

#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

/** 
 * @brief Class GeantFactoryStore
 * @details  Store of factories holding one factory per user type.
 */
class GeantFactoryStore : public TObject {
private:
  Int_t fNclients;                      /** Number of thread clients */
  Int_t fNFactories;                    /** Number of factories stored */
  Int_t fCapacity;                      /** Store capacity */
  const void **fTypes;                  /** Factory types */
  void **fFactories;                    /** Stored factories */
  static GeantFactoryStore *fgInstance; /** Static instance of the store */
  
  /**
   * @brief Constructor GeantFactoryStore
   * 
   * @param nclients Number of thread clients 
   */
  GeantFactoryStore(Int_t nclients);
  
  /**
   * @brief Function of removal factory
   * 
   * @param islot Slot ID
   */
  void RemoveFactory(Int_t islot);

  /**
   * @brief Constructor GeantFactoryStore
   * 
   * @todo Still needs to be implemented
   */
  GeantFactoryStore(const GeantFactoryStore &);

  /**
   * @brief Implementation of operator=
   * 
   * @todo Still needs to be implemented
   */
  GeantFactoryStore &operator=(const GeantFactoryStore &);
public:
  
  /** @brief Destructor GeantFactoryStore */
  virtual ~GeantFactoryStore();
  
  /**
   * @brief  Function that creates object of instance for factory
   * 
   * @param nclients Number of thread clients (by default 1)
   */
  static GeantFactoryStore *Instance(Int_t nclients = 1);
  
  /**
   * @brief Templated function that return factory object
   * 
   * @param blocksize Size of block
   */
  template <class T> GeantFactory<T> *GetFactory(Int_t blocksize);
  
  /** @brief Function that provides deletion of factory */
  template <class T> void DeleteFactory();
};

/**
 * @details Returns existing factory for the user type or create new one.
 * @return Factory object
 */
template <class T> GeantFactory<T> *GeantFactoryStore::GetFactory(Int_t blocksize) {
  const std::type_info *type = &typeid(T);
  for (Int_t i = 0; i < fNFactories; i++) {
    if ((const std::type_info *)fTypes[i] == type)
      return (GeantFactory<T> *)fFactories[i];
  }
  GeantFactory<T> *factory = new GeantFactory<T>(fNclients, blocksize);
  TThread::Lock();
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
  TThread::UnLock();
  return factory;
}

/**
 * @todo Add some small details for doxygen
 */
template <class T> void GeantFactoryStore::DeleteFactory() {
  const std::type_info *type = &typeid(T);
  for (Int_t i = 0; i < fNFactories; i++) {
    if ((const std::type_info *)fTypes[i] == type) {
      GeantFactory<T> *factory = (GeantFactory<T> *)fFactories[i];
      delete factory;
      RemoveFactory(i);
      return;
    }
  }
}
#endif
