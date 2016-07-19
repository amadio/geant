#ifndef TRANSPORTTASK
#define TRANSPORTTASK

#include "tbb/task.h"

class TransportTask : public tbb::task
{
private:
  int fNbaskets;

public:
  TransportTask (int nbaskets);
  ~TransportTask ();

  tbb::task* execute ();

};

#endif // TRANSPORTTASK
