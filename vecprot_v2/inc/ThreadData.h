#ifndef THREADDATA
#define THREADDATA

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif
#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif
#include "TThread.h"
#include "GeantFactoryStore.h"
#include "TThreadMergingFile.h"

class ThreadData {
protected:
  static ThreadData *fgInstance;
public:
  // Data per thread
  int fNthreads;
  Geant::TThreadMergingFile **fFiles;
  TTree **fTrees;
  GeantBlock<MyHit> **fData;
  GeantBasketMgr **fPrioritizers;
  GeantFactory<MyHit> **fMyhitFactories;

  ThreadData (int nthreads);
  ~ThreadData ();

  static ThreadData *Instance(int nthreads);

};

#endif //THREADDATA
