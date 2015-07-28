#ifndef ROOT_TNudyDB
#define ROOT_TNudyDB

#include <THashTable.h>
#include <TFile.h>
#include <TList.h>

#include "TNudyCore.h"
#include "TNudyLibrary.h"


class TNudyDB : public TNamed {
 public:
  //   TNudyDB(const char *name, const char *title, const char *file);
  TNudyDB(const Char_t *name, const Char_t *title, const char *file);// DEBUG
  virtual ~TNudyDB();
  void AddLibrary(const Char_t *name, const char* file); //Add a Library to a Nudy Database from a RENDF file
  void RemoveLibrary(const Char_t *name); //Remove a Library from a Nudy Databse
  TList* GetEntries();  //Get a List of Libraries in the Nudy Database
  TFile* GetDBFile() {return fDB;} //Get a pointer to the database file
 private:
  
  TFile *fDB; //! Database data file
  
  ClassDef(TNudyDB,1)
};
#endif
