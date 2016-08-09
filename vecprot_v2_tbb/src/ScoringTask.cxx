#include "Geant/Error.h"

#include "ScoringTask.h"

#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif


ScoringTask::ScoringTask (Geant::GeantTaskData *td): fTd(td) { }

ScoringTask::~ScoringTask () { }

tbb::task* ScoringTask::execute ()
{
  return NULL;
}
