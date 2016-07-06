#ifndef TRANSPORTTASK
#define TRANSPORTTASK

#include "tbb/task.h"

class TransportTask : public tbb::task
{
private:

public:
  TransportTask ();
  ~TransportTask ();

  tbb::task* execute ();

};

#endif // TRANSPORTTASK
