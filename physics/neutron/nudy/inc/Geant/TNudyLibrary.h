#ifndef ROOT_TNudyLibrary
#define ROOT_TNudyLibrary

#include "TNamed.h"
#include "THashTable.h"
#include "TParticlePDG.h"
class TNudyEndfTape;
class TParticlePDG;
class TNudySubLibrary;
class TGeoElementRN;

class TNudyLibrary : public TNamed {
public:
  TNudyLibrary();                                     // Default Constructor
  TNudyLibrary(const char *name, const char *title);  // Constructor
  void ReadTape(TNudyEndfTape *tape);                 // Read data from tape into current file
  virtual ~TNudyLibrary();                            // Destructor
  TNudySubLibrary *AddSubLib(TParticlePDG *particle); // Add a sublibrary for particle
  THashTable *GetSubLibs() { return fSubLib; }        // Get all sublibraries
  TNudySubLibrary *GetSubLib(TParticlePDG *particle)
  {
    return (TNudySubLibrary *)fSubLib->FindObject(particle->GetName());
  } // Get sublibrary by particle
  TNudySubLibrary *SetSubLib(TParticlePDG *particle) { return fCurSubLib = GetSubLib(particle); }
  TNudySubLibrary *GetSubLib() { return fCurSubLib; }
  void ListModels();
  bool IsHandled(TParticlePDG *particle, TGeoElementRN *target, unsigned long temp);

private:
  THashTable *fSubLib;         // Sub-Libaries storing the MAT-MT information for each particle
  TNudySubLibrary *fCurSubLib; //! Current Sublibrary

#ifdef USE_ROOT
  ClassDef(TNudyLibrary, 1)
#endif
};
#endif
