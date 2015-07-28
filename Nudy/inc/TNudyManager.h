#ifndef ROOT_TNudyManager
#define ROOT_TNudyManager

#include "TNudyCore.h"
#include "TNudyDB.h"
#include "TNudyLibrary.h"
#include "TNudySubLibrary.h"
#include "TVNudyModel.h"
#include "TNudyENDF.h"

class TNudyManager : public TNamed {
 protected:
  static TNudyManager* fgInstance;
  TNudyManager();
  TNudyManager(const TNudyManager& manager): TNamed(manager) {}
  TNudyManager& operator=(const TNudyManager& manager)
    {if(this!=&manager){TNamed::operator=(manager);}return *this;}
  THashTable *fNudyDB; 
  THashTable *fLibrary;
  TNudyDB *fCurNudyDB;
  TNudyLibrary *fCurLibrary;
  TNudyCore *fCore;
  TBtree *fResult;

 public:
  ~TNudyManager();//Public Destructor
  static TNudyManager* Instance();//Returns Instance of TNudyManager
  TNudyLibrary* LoadLibrary(const char* memLibName, const char* diskLibName,const char* sublib = NULL, TGeoElementRN *mat = NULL, Reaction_t reac = kNoReaction,ULong_t temp = 0);
  THashTable * GetDatabases() {return fNudyDB;}
  THashTable * GetLibraries() {return fLibrary;}
  TNudyDB* OpenDatabase(const char* name, const char* file);
  TNudyDB* SetDatabase(const char* name);
  TNudyDB* SetDatabase(const TNudyDB* db);
  TNudyDB* GetDatabase(){return fCurNudyDB;}
  TNudyLibrary* GetLibrary(const char *name){return (TNudyLibrary*)fLibrary->FindObject(name);}
  TNudyLibrary* SetLibrary(TNudyLibrary* lib){return (fCurLibrary = (TNudyLibrary*)fLibrary->FindObject(lib));}
  TNudyLibrary* SetLibrary(const char *name){return (fCurLibrary = GetLibrary(name));}
  void DumpTape(const char* rendf, Int_t debug);
  void SetSubLib(TParticlePDG *proj){fCurLibrary->SetSubLib(proj);}
  void SetSubLib(const char *name){fCurLibrary->SetSubLib(TNudyCore::Instance()->GetParticlePDG(name));}
  void ProcessTape(const char* endf, const char* rendf);
  Int_t CloseDatabase(const char* name = NULL);
  void AddEndfLibrary(const char* name, const char *endf);
  void AddLibrary(const char *name, const char *file){if(fCurNudyDB) fCurNudyDB->AddLibrary(name,file);}
  void RemoveLibrary(const char *name){if(fCurNudyDB) fCurNudyDB->RemoveLibrary(name);}
  void ListModels();
  TVNudyModel* GetModel(const Int_t a, const Int_t z , const Int_t iso, const Int_t reaction,const ULong_t temp, const char* particleName = NULL);
  TBtree *GetAllModels(const Int_t a = 0,const Int_t z = 0 ,const Int_t iso = 0,const Int_t reaction = kNoReaction,const ULong_t temp = 0,const char *particleName = NULL);

  ClassDef(TNudyManager,1)

};
#endif 
