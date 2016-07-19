#ifndef FEEDERTASK
#define FEEDERTASK

#include "tbb/task.h"

class FeederTask : public tbb::task
{
private:
  int *fNbaskets;

public:
  FeederTask (int *nbaskets);
  ~FeederTask ();

  tbb::task* execute ();

};

#endif //FEEDERTASK
