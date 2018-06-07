#include "Geant/TNudyDB.h"
#include "Geant/TNudyEndfTape.h"
#include "TFile.h"
#include "Geant/TNudyLibrary.h"

using namespace Nudy;

#ifdef USE_ROOT
ClassImp(TNudyDB)
#endif

//______________________________________________________________________________
TNudyDB::TNudyDB(const char *name, const char *title, const char *file)
: fName(name), fTitle(title)
{
  // Constructor for TNudyDB, opens database storage file
  // Try to open file for update
  fDB = new TFile(file, "UPDATE");
  if (fDB->IsOpen()) {
    // if file exists open for reading
    printf("Database Opened\n");
    // Create new library hash in memory
  } else {
    printf("Creating new Database\n");
    fDB->Close();
    fDB = new TFile(file, "RECREATE");
    // **Check if creation is successful
  }
  //  fDB->Close();
}

//______________________________________________________________________________
TNudyDB::~TNudyDB()
{
  // Destructor for TNudyDB, Deletes all libraries in memory
  
  // Close open database file
  if (fDB->IsOpen()) {
    fDB->Close();
  }
  // Delete database file in memory
  delete fDB;
}

//______________________________________________________________________________
void TNudyDB::AddLibrary(const char *name, const char *file)
{
  // Open endf root file and store tape
  // Check if file exists
  TFile *endf = new TFile(file, "OLD");
  if (endf->IsZombie()) {
    //    TNudyManager::Instance()->Log("Error opening ENDF file");
    Error("AddLibrary", "Cannot open file %s", file);
    return;
  }
  Nudy::TNudyEndfTape *newTape = new Nudy::TNudyEndfTape();
  // Check if tape is present
  if (newTape) {
    newTape->Read(endf->GetListOfKeys()->First()->GetName());
    if (!newTape) {
      Error("AddLibrary", "Cannot find ENDF Tape in file %s", file);
      return;
    }
  } else {
    //    TNudyManager::Instance()->Log("Error opening ENDF Tape");
    Error("AddLibrary", "Could not create new ENDF Tape");
    return;
  }
  // Create tempory library in memory to process and store data
  Nudy::TNudyLibrary *lib = new Nudy::TNudyLibrary(name, file);
  
  // Set current file to Database file
  gFile = fDB;
  gDirectory->cd("/");
  if (!fDB->GetDirectory(name)) {
    // If library is not stored in database create new library
    printf("Creating new Library\n");
    fDB->mkdir(name);
  } else {
    // If library exists open it for updating
    printf("Updating Library %s\n", name);
    //    lib = (Nudy::TNudyLibrary*)fDB->Get(name);
    // lib->Print();
  }
  
  // Go to Library folder
  fDB->cd(name);
  // gDirectory->cd(name);
  //  printf("%s",fDB->pwd());
  
  // Read and Process tape
  lib->ReadTape(newTape);
  
  // Close file
  endf->Close();
  
  // Return to root directory of database file
  gDirectory->cd("/");
  fDB->cd("/");
  // Delete endf file , temporary library and tape from memory
  delete endf;
  delete lib;
  printf("Deleting Tape");
  delete newTape;
}

//______________________________________________________________________________
TList *TNudyDB::GetEntries()
{
  return fDB->GetListOfKeys();
}

//______________________________________________________________________________
void TNudyDB::RemoveLibrary(const char *name)
{
  fDB->Delete(Form("%s;*", name));
}
