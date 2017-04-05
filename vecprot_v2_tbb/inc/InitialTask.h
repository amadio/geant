#ifndef GEANT_TBB_INITIAL_TASK
#define GEANT_TBB_INITIAL_TASK

#include "tbb/task.h"
#include "GeantPropagator.h"

class InitialTask : public tbb::task
{
public:

private:
  Geant::GeantPropagator *fPropagator;

public:
  InitialTask(Geant::GeantPropagator *prop) : fPropagator(prop) {}
  ~InitialTask() {}

  tbb::task* execute();

};

#endif //GEANT_TBB_INITIAL_TASK
