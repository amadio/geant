#ifndef GEANT_TBB_FLOW_CONTROLLER_TASK
#define GEANT_TBB_FLOW_CONTROLLER_TASK

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

#endif //GEANT_TBB_FLOW_CONTROLLER_TASK
