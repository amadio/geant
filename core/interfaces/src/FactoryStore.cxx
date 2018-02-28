#include "Geant/FactoryStore.h"

FactoryStore *FactoryStore::fgInstance = 0;

/**
 * @details Creates a factory store that can serve a
 * given number of thread clients.
 */
//______________________________________________________________________________
FactoryStore::FactoryStore(int nclients)
    : fNclients(nclients), fNFactories(0), fCapacity(100), fTypes(0), fFactories(0), fMutex()
{
  // Private constructor
  fTypes     = new const void *[fCapacity];
  fFactories = new void *[fCapacity];
  memset(fTypes, 0, fCapacity * sizeof(void *));
  memset(fFactories, 0, fCapacity * sizeof(void *));
}

/**
 * @details Returns the pointer to the factory store singleton. If not
 * existing create one.
 */
//______________________________________________________________________________
FactoryStore *FactoryStore::Instance(int nclients)
{
  // Returns singleton for the factory store.
  if (fgInstance) return fgInstance;
  fgInstance = new FactoryStore(nclients);
  return fgInstance;
}

/**
 * @Delete all stored types and factories, ressetting the singleton pointer
 */
//______________________________________________________________________________
FactoryStore::~FactoryStore()
{
  // Destructor
  delete[] fTypes;
  //   for (int i=0; i<fNFactories; i++) delete fFactories[i];
  delete[] fFactories;
  fgInstance = 0;
}

/**
 * @Delete a single factory at the given slot.
 */
//______________________________________________________________________________
void FactoryStore::RemoveFactory(int islot)
{
  // Remove a factory from the store. Called by the destructor of the factory.
  fNFactories--;
  if (islot == fNFactories) return;
  memmove(&fFactories[islot], &fFactories[islot + 1], (fNFactories - islot) * sizeof(void *));
  memmove(&fTypes[islot], &fTypes[islot + 1], (fNFactories - islot) * sizeof(void *));
}
