#ifndef ROOT_TNudyManager
#define ROOT_TNudyManager

#include "TNamed.h"
#include "TBtree.h"
class TVNudyModel;
#include "THashTable.h"
#include "TNudyDB.h"
#include "TNudyLibrary.h"
#include "TNudyCore.h"
class TGeoElementRN;
class TParticlePDG;
#include "TNudyTypes.h"

class TNudyManager : public TNamed {
protected:
  static TNudyManager *fgInstance;
  TNudyManager();
  TNudyManager(const TNudyManager &manager) : TNamed(manager) {}
  TNudyManager &operator=(const TNudyManager &manager) {
    if (this != &manager) {
      TNamed::operator=(manager);
    }
    return *this;
  }
  THashTable *fNudyDB;
  THashTable *fLibrary;
  TNudyDB *fCurNudyDB;
  TNudyLibrary *fCurLibrary;
  TNudyCore *fCore;
  TBtree *fResult;

public:
  ~TNudyManager();                 // Public Destructor
  static TNudyManager *Instance(); // Returns Instance of TNudyManager
  TNudyLibrary *LoadLibrary(const char *memLibName, const char *diskLibName, const char *sublib = NULL,
                            TGeoElementRN *mat = NULL, Reaction_t reac = kNoReaction, unsigned long temp = 0);
  THashTable *GetDatabases() { return fNudyDB; }
  THashTable *GetLibraries() { return fLibrary; }
  TNudyDB *OpenDatabase(const char *name, const char *file);
  TNudyDB *SetDatabase(const char *name);
  TNudyDB *SetDatabase(const TNudyDB *db);
  TNudyDB *GetDatabase() { return fCurNudyDB; }
  TNudyLibrary *GetLibrary(const char *name) { return (TNudyLibrary *)fLibrary->FindObject(name); }
  TNudyLibrary *SetLibrary(TNudyLibrary *lib) { return (fCurLibrary = (TNudyLibrary *)fLibrary->FindObject(lib)); }
  TNudyLibrary *SetLibrary(const char *name) { return (fCurLibrary = GetLibrary(name)); }
  void DumpTape(const char *rendf, int debug);
  void SetSubLib(TParticlePDG *proj) { fCurLibrary->SetSubLib(proj); }
  void SetSubLib(const char *name) { fCurLibrary->SetSubLib(TNudyCore::Instance()->GetParticlePDG(name)); }
  void ProcessTape(const char *endf, const char *rendf);
  int CloseDatabase(const char *name = NULL);
  void AddEndfLibrary(const char *name, const char *endf);
  void AddLibrary(const char *name, const char *file) {
    if (fCurNudyDB)
      fCurNudyDB->AddLibrary(name, file);
  }
  void RemoveLibrary(const char *name) {
    if (fCurNudyDB)
      fCurNudyDB->RemoveLibrary(name);
  }
  void ListModels();
  TVNudyModel *GetModel(const int a, const int z, const int iso, const int reaction, const unsigned long temp,
                        const char *particleName = NULL);
  TBtree *GetAllModels(const int a = 0, const int z = 0, const int iso = 0, const int reaction = kNoReaction,
                       const unsigned long temp = 0, const char *particleName = NULL);

#ifdef USE_ROOT
  ClassDef(TNudyManager, 1)
#endif
};
#endif
