#ifndef FLOWCONTROLLERTASK
#define FLOWCONTROLLERTASK

#include "GeantTaskData.h"
#include "ThreadData.h"
#include "FeederTask.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "TThreadMergingFile.h"

#ifdef GEANT_TBB
#include "tbb/task.h"
#endif

class FlowControllerTask : public tbb::task
{
private:
  Geant::GeantTaskData *fTd;
  bool fStarting;

public:
  FlowControllerTask (Geant::GeantTaskData *td, bool starting);
  ~FlowControllerTask ();

  tbb::task* execute ();

};

#endif //FLOWCONTROLLERTASK
