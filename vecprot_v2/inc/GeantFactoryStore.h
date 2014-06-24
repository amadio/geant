#ifndef GEANT_FACTORYSTORE
#define GEANT_FACTORYSTORE

#include <Rtypes.h>
#include <typeinfo>

#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

//______________________________________________________________________________
// GeantFactoryStore - Store of factories holding one factory per user type.
//______________________________________________________________________________

class GeantFactoryStore : public TObject {
private:
   Int_t                fNclients;        // Number of thread clients
   Int_t                fNFactories;      // Number of factories stored
   Int_t                fCapacity;        // Store capacity
   void               **fTypes;           // Factory types
   void               **fFactories;       // Stored factories
   static GeantFactoryStore *fgInstance;  // Static instance of the store
   GeantFactoryStore(Int_t nclients);
   void                 RemoveFactory(Int_t islot);
   GeantFactoryStore(const GeantFactoryStore&); // Not implemented
   GeantFactoryStore &operator=(const GeantFactoryStore&); // Not implemented
public:
   virtual ~GeantFactoryStore();
   
   static GeantFactoryStore *
                        Instance(Int_t nclients=1);
   template <class T>
   GeantFactory<T>     *GetFactory(Int_t blocksize);
   template <class T>
   void                 DeleteFactory();
};

template <class T>
GeantFactory<T> *GeantFactoryStore::GetFactory(Int_t blocksize)
{
// Returns existing factory for the user type or create new one.
   const std::type_info *type = &typeid(T);
   for (Int_t i=0; i<fNFactories; i++) {
      if ((const std::type_info*)fTypes[i]==type) return (GeantFactory<T>*)fFactories[i];
   }
   GeantFactory<T> *factory = new GeantFactory<T>(fNclients, blocksize);
   TThread::Lock();
   if (fNFactories==fCapacity) {
      // Resize arrays
      void **types = new void*[2*fCapacity];
      memset(types, 0, 2*fCapacity*sizeof(void*));
      memcpy(types, fTypes, fNFactories*sizeof(void*));
      void **factories = new void*[2*fCapacity];
      memset(factories, 0, 2*fCapacity*sizeof(void*));
      memcpy(factories, fFactories, fNFactories*sizeof(void*));
      delete [] fTypes;
      fTypes = types;
      delete [] fFactories;
      fFactories = factories;
      fCapacity *= 2;
   }   
   fTypes[fNFactories] = (void*)type;
   fFactories[fNFactories++] = factory;
   TThread::UnLock();
   return factory;
}

template <class T>
void GeantFactoryStore::DeleteFactory()
{
// Deletes the factory for this type
   const std::type_info *type = &typeid(T);
   for (Int_t i=0; i<fNFactories; i++) {
      if ((const std::type_info*)fTypes[i]==type) {
         GeantFactory<T> *factory = (GeantFactory<T>*)fFactories[i];
         delete factory;
         RemoveFactory(i);
         return;
      }
   }
}
#endif
