#ifndef FLOWCONTROLLERTASK
#define FLOWCONTROLLERTASK

#include "GeantTaskData.h"
#include "tbb/task.h"

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
