#ifndef TNudyDB_H
#define TNudyDB_H

#include <string>
#include "TObject.h"
class TFile;

using std::string;

namespace Nudy {
class TNudyDB : public TObject {
public:
  TNudyDB(const char *name, const char *title, const char *file); // DEBUG
  virtual ~TNudyDB();
  void AddLibrary(const char *name, const char *file); // Add a Library to a Nudy Database from a RENDF file
  void RemoveLibrary(const char *name);                // Remove a Library from a Nudy Databse
  TList *GetEntries();                                 // Get a List of Libraries in the Nudy Database
  TFile *GetDBFile() { return fDB; }                   // Get a pointer to the database file
private:
  string fName;
  string fTitle;
  TFile *fDB; //! Database data file

#ifdef USE_ROOT
  ClassDef(TNudyDB, 1)
#endif
};

} // namespace
#endif
