#include "GeantFactoryStore.h"

GeantFactoryStore *GeantFactoryStore::fgInstance = 0;

//______________________________________________________________________________
GeantFactoryStore::GeantFactoryStore(Int_t nclients)
                  :fNclients(nclients),
                   fNFactories(0),
                   fCapacity(100),
                   fTypes(0),
                   fFactories(0)
{
// Private constructor
   fTypes = new void*[fCapacity];
   fFactories = new void*[fCapacity];
   memset(fTypes, 0, fCapacity*sizeof(void*));
   memset(fFactories, 0, fCapacity*sizeof(void*));
}

//______________________________________________________________________________
GeantFactoryStore *GeantFactoryStore::Instance(Int_t nclients)
{
// Returns singleton for the factory store.
   if (fgInstance) return fgInstance;
   fgInstance = new GeantFactoryStore(nclients);
   return fgInstance;
}

//______________________________________________________________________________
GeantFactoryStore::~GeantFactoryStore()
{
// Destructor
   delete [] fTypes;
//   for (Int_t i=0; i<fNFactories; i++) delete fFactories[i];
   delete [] fFactories;
   fgInstance = 0;
}

//______________________________________________________________________________
void GeantFactoryStore::RemoveFactory(Int_t islot)
{
// Remove a factory from the store. Called by the destructor of the factory.
   fNFactories--;
   if (islot==fNFactories) return;
   memmove(&fFactories[islot], &fFactories[islot+1], (fNFactories-islot)*sizeof(void*));
   memmove(&fTypes[islot], &fTypes[islot+1], (fNFactories-islot)*sizeof(void*));
}
