#ifndef FLOWCONTROLLERTASK
#define FLOWCONTROLLERTASK

#include "tbb/task.h"

class FlowControllerTask : public tbb::task
{
private:
  bool fStarting;
public:
  FlowControllerTask (bool starting);
  ~FlowControllerTask ();

  tbb::task* execute ();

};

#endif //FLOWCONTROLLERTASK
