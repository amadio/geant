#ifndef ROOT_TNudyLibrary
#define ROOT_TNudyLibrary

#include "TNamed.h"
#include "THashTable.h"
#include "TParticlePDG.h"

namespace Nudy {
class TNudyEndfTape;
class TNudySubLibrary;
}

class TParticlePDG;
class TGeoElementRN;

namespace Nudy {
class TNudyLibrary : public TNamed {
public:
  TNudyLibrary();                                           // Default Constructor
  TNudyLibrary(const char *name, const char *title);        // Constructor
  void ReadTape(Nudy::TNudyEndfTape *tape);                 // Read data from tape into current file
  virtual ~TNudyLibrary();                                  // Destructor
  Nudy::TNudySubLibrary *AddSubLib(TParticlePDG *particle); // Add a sublibrary for particle
  THashTable *GetSubLibs() { return fSubLib; }              // Get all sublibraries
  Nudy::TNudySubLibrary *GetSubLib(TParticlePDG *particle)
  {
    return (Nudy::TNudySubLibrary *)fSubLib->FindObject(particle->GetName());
  } // Get sublibrary by particle
  Nudy::TNudySubLibrary *SetSubLib(TParticlePDG *particle) { return fCurSubLib = GetSubLib(particle); }
  Nudy::TNudySubLibrary *GetSubLib() { return fCurSubLib; }
  void ListModels();
  bool IsHandled(TParticlePDG *particle, TGeoElementRN *target, unsigned long temp);

private:
  THashTable *fSubLib;               // Sub-Libaries storing the MAT-MT information for each particle
  Nudy::TNudySubLibrary *fCurSubLib; //! Current Sublibrary

#ifdef USE_ROOT
  ClassDef(TNudyLibrary, 1)
#endif
};

} // namespace
#endif
