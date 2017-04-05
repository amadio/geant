#ifndef GEANT_TBB_FLOW_CONTROLLER_TASK
#define GEANT_TBB_FLOW_CONTROLLER_TASK

#include "GeantTaskData.h"
#include "tbb/task.h"

class FlowControllerTask : public tbb::task
{
private:
  Geant::GeantTaskData *fTd;
  bool fStarting;
  bool fForcedStop;

public:
  FlowControllerTask (Geant::GeantTaskData *td, bool starting, bool forcedStop=false);
  ~FlowControllerTask ();

  tbb::task *SpawnFeederTask();
  tbb::task *SpawnTransportTask();

  tbb::task *execute();

};

#endif //GEANT_TBB_FLOW_CONTROLLER_TASK
