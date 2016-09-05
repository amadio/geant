#ifndef GEANT_TBB_INITIAL_TASK
#define GEANT_TBB_INITIAL_TASK

#include "tbb/task.h"
#include "GeantTaskData.h"

class GeantPropagator;

class InitialTask : public tbb::task
{
private:
  GeantPropagator *fPropagator;

public:
  InitialTask(GeantPropagator *prop) : fPropagator(prop) {}
  ~InitialTask() {}

  tbb::task* execute();

};

#endif //GEANT_TBB_INITIAL_TASK
