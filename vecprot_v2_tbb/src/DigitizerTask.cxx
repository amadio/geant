#include "Geant/Error.h"

#include "DigitizerTask.h"

#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif


DigitizerTask::DigitizerTask (Geant::GeantTaskData *td): fTd(td) { }

DigitizerTask::~DigitizerTask () { }

tbb::task* DigitizerTask::execute ()
{
  return NULL;
}
