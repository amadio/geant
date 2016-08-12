#ifndef INITIALTASK
#define INITIALTASK

#include "tbb/task.h"

class InitialTask : public tbb::task
{
private:

public:
  InitialTask() {}
  ~InitialTask() {}

  tbb::task* execute();

};

#endif //INITIALTASK
