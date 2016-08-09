#include "Geant/Error.h"

#include "IOTask.h"

#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif


IOTask::IOTask (Geant::GeantTaskData *td): fTd(td) { }

IOTask::~IOTask () { }

tbb::task* IOTask::execute ()
{
  return NULL;
}
