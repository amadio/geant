#include "ThreadData.h"
#include "Geant/Error.h"

using namespace Geant;

ThreadData *ThreadData::fgInstance = 0;

ThreadData::ThreadData(int nthreads): fNthreads(nthreads){
  // initialize data per thread
  fFiles = new TThreadMergingFile *[fNthreads];
  fTrees = new TTree *[fNthreads];
  fData = new GeantBlock<MyHit> *[fNthreads];
  fPrioritizers = new GeantBasketMgr *[fNthreads];
  fMyhitFactories = new GeantFactory<MyHit> *[fNthreads];
}

ThreadData::~ThreadData() {
  // Destructor.
  delete[] fFiles;
  delete[] fTrees;
  delete[] fData;
  delete[] fPrioritizers;
  delete[] fMyhitFactories;
}

ThreadData *ThreadData::Instance(int nthreads) {
  // Return singleton instance.
  if (fgInstance)
    return fgInstance;
  if (!nthreads) {
    ::Error("ThreadData::Instance", "No instance yet so you should provide number of threads.");
    return 0;
  }
  return new ThreadData(nthreads);
}
