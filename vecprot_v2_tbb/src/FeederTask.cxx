#include "FeederTask.h"

#include "GeantConfig.h"
#include "GeantRunManager.h"
#include "GeantScheduler.h"

#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif


FeederTask::FeederTask (Geant::GeantTaskData *td, bool starting): fTd(td), fStarting(starting) { }

FeederTask::~FeederTask () { }

tbb::task* FeederTask::execute ()
{
  GeantRunManager *runmgr = fTd->fPropagator->fRunMgr;
  int nbaskets = runmgr->Feeder(fTd);
  if (fStarting && (nbaskets == 0)) fTd->fPropagator->fWMgr->GetScheduler()->GarbageCollect(fTd, true);
  tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
  TransportTask & transportTask = *new(cont.allocate_child()) TransportTask( fTd );
  return & transportTask;
}
