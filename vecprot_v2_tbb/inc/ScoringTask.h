#ifndef GEANT_TBB_SCORING_TASK
#define GEANT_TBB_SCORING_TASK

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

#endif //GEANT_TBB_SCORING_TASK
