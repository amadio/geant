#include "ThreadData.h"
#include "Geant/Error.h"

#ifdef USE_ROOT
#include "TTree.h"
#include "TThreadMergingFile.h"
#endif

ThreadData *ThreadData::fgInstance = 0;

ThreadData::ThreadData(int nthreads): fNthreads(nthreads){
  // initialize data per thread
#ifdef USE_ROOT
  fFiles = new Geant::TThreadMergingFile *[fNthreads];
  fTrees = new TTree *[fNthreads];
#endif
  fData = new GeantBlock<MyHit> *[fNthreads];
  fMyhitFactories = new GeantFactory<MyHit> *[fNthreads];

  // instance
  fgInstance = this;
}

ThreadData::~ThreadData() {
  // Destructor.
#ifdef USE_ROOT
  delete[] fFiles;
  delete[] fTrees;
#endif
  delete[] fData;
  delete[] fMyhitFactories;
}

ThreadData *ThreadData::Instance(int nthreads) {
  // Return singleton instance.
  if (fgInstance)
    return fgInstance;
  if (!nthreads) {
    Geant::Error("ThreadData::Instance", "No instance yet so you should provide number of threads.");
    return 0;
  }
  return new ThreadData(nthreads);
}
