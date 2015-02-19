#include "GeantFactoryStore.h"

GeantFactoryStore *GeantFactoryStore::fgInstance = 0;

/**
 * @details Creates a factory store that can serve a
 * given number of thread clients.
 */
//______________________________________________________________________________
GeantFactoryStore::GeantFactoryStore(int nclients)
    : fNclients(nclients), fNFactories(0), fCapacity(100), fTypes(0), fFactories(0), fMutex() {
  // Private constructor
  fTypes = new const void *[fCapacity];
  fFactories = new void *[fCapacity];
  memset(fTypes, 0, fCapacity * sizeof(void *));
  memset(fFactories, 0, fCapacity * sizeof(void *));
}

/**
 * @details Returns the pointer to the factory store singleton. If not
 * existing create one.
 */
//______________________________________________________________________________
GeantFactoryStore *GeantFactoryStore::Instance(int nclients) {
  // Returns singleton for the factory store.
  if (fgInstance)
    return fgInstance;
  fgInstance = new GeantFactoryStore(nclients);
  return fgInstance;
}

/**
 * @Delete all stored types and factories, ressetting the singleton pointer
 */
//______________________________________________________________________________
GeantFactoryStore::~GeantFactoryStore() {
  // Destructor
  delete [] fTypes;
  //   for (Int_t i=0; i<fNFactories; i++) delete fFactories[i];
  delete [] fFactories;
  fgInstance = 0;
}

/**
 * @Delete a single factory at the given slot.
 */
//______________________________________________________________________________
void GeantFactoryStore::RemoveFactory(int islot) {
  // Remove a factory from the store. Called by the destructor of the factory.
  fNFactories--;
  if (islot == fNFactories)
    return;
  memmove(&fFactories[islot], &fFactories[islot + 1], (fNFactories - islot) * sizeof(void *));
  memmove(&fTypes[islot], &fTypes[islot + 1], (fNFactories - islot) * sizeof(void *));
}
