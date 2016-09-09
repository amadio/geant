#ifndef GEANT_TBB_THREAD_DATA
#define GEANT_TBB_THREAD_DATA

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif
#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif
#include "TThread.h"
#include "GeantFactoryStore.h"
#include "TThreadMergingFile.h"

#ifdef USE_ROOT
class TTree;
class TThreadMergingFile;
#endif
class GeantBasketMgr;
template <typename T>
class GeantFactory;

class ThreadData {
protected:
  static ThreadData *fgInstance;
public:
  // Data per thread
  int fNthreads;
#ifdef USE_ROOT
  Geant::TThreadMergingFile **fFiles;
  TTree **fTrees;
#endif
  GeantBlock<MyHit> **fData;
  GeantFactory<MyHit> **fMyhitFactories;

  ThreadData (int nthreads);
  ~ThreadData ();

  static ThreadData *Instance(int nthreads);

};

#endif //GEANT_TBB_THREAD_DATA
