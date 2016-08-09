#ifndef SCORINGTASK
#define SCORINGTASK

#include "GeantTaskData.h"
#include "tbb/task.h"

class ScoringTask : public tbb::task
{
private:
  Geant::GeantTaskData *fTd;

public:
  ScoringTask (Geant::GeantTaskData *td);
  ~ScoringTask ();

  tbb::task* execute ();

};

#endif //SCORINGTASK
