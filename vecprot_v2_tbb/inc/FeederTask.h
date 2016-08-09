#ifndef FEEDERTASK
#define FEEDERTASK

#include "GeantTaskData.h"
#include "tbb/task.h"

class FeederTask : public tbb::task
{
private:
  Geant::GeantTaskData *fTd;
  int *fNbaskets;

public:
  FeederTask (Geant::GeantTaskData *td, int *nbaskets);
  ~FeederTask ();

  tbb::task* execute ();

};

#endif //FEEDERTASK
