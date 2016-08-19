#ifndef GEANT_TBB_INITIAL_TASK
#define GEANT_TBB_INITIAL_TASK

#include "tbb/task.h"

class InitialTask : public tbb::task
{
private:

public:
  InitialTask() {}
  ~InitialTask() {}

  tbb::task* execute();

};

#endif //GEANT_TBB_INITIAL_TASK
