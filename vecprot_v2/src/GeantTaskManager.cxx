#include "GeantTaskManager.h"
#include "TThread.h"

ClassImp(GeantTaskManager)

    //______________________________________________________________________________
    GeantTaskManager::GeantTaskManager(const char *name)
    : TNamed(name, ""), fNthreads(0), fNactive(0), fQueue() {
  // Constructor.
}

//______________________________________________________________________________
GeantTaskManager::~GeantTaskManager() {
  // Destructor
  StopWorkers();
}

//______________________________________________________________________________
Int_t GeantTaskManager::AddWorkers(Int_t nworkers) {
  // Add new threads to the workers pool.
  for (Int_t i = 0; i < nworkers; i++) {
    TThread *t = new TThread(GeantTaskManager::ProcessLoop);
    fThreads->push_back(t);
    t->Run();
  }
}
