#ifndef FEEDERTASK
#define FEEDERTASK

#include <atomic>

#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "GeantEvent.h"

#include "tbb/task.h"

using namespace tbb;

class FeederTask : public task
{
private:
  Geant::GeantTaskData *fTd;
  int * fNbaskets;

public:
  FeederTask (Geant::GeantTaskData *td, int *nbaskets);
  ~FeederTask ();

  task* execute ();

};

#endif //FEEDERTASK
