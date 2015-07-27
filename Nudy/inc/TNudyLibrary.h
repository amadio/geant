#ifndef ROOT_TNudyLibrary
#define ROOT_TNudyLibrary

#include <THashList.h>
#include <THashTable.h>
#include "TNudyCore.h"
#include "TNudySubLibrary.h"
#include "TNudyEndfMat.h"
#include "TNudyEndfTape.h"

class TNudyLibrary : public TNamed {
public:
  TNudyLibrary(); // Default Constructor
  TNudyLibrary(const char *name, const char *title); // Constructor
  void ReadTape(TNudyEndfTape *tape); // Read data from tape into current file
  virtual ~TNudyLibrary(); // Destructor
  TNudySubLibrary *AddSubLib(TParticlePDG *particle); // Add a sublibrary for particle
  THashTable *GetSubLibs() { return fSubLib; } // Get all sublibraries
  TNudySubLibrary *GetSubLib(TParticlePDG *particle) {
    return (TNudySubLibrary *)fSubLib->FindObject(particle->GetName());
  } // Get sublibrary by particle
  TNudySubLibrary *SetSubLib(TParticlePDG *particle) { return fCurSubLib = GetSubLib(particle); }
  TNudySubLibrary *GetSubLib() { return fCurSubLib; }
  void ListModels();
  Bool_t IsHandled(TParticlePDG *particle, TGeoElementRN *target, ULong_t temp);

private:
  THashTable *fSubLib;         // Sub-Libaries storing the MAT-MT information for each particle
  TNudySubLibrary *fCurSubLib; //! Current Sublibrary

  ClassDef(TNudyLibrary, 1)
};
#endif
