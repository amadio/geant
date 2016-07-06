#ifndef FLOWCONTROLLERTASK
#define FLOWCONTROLLERTASK

#include "tbb/task.h"

class FlowControllerTask : public tbb::task
{
private:

public:
  FlowControllerTask ();
  ~FlowControllerTask ();

  tbb::task* execute ();

};

#endif //FLOWCONTROLLERTASK
