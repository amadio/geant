#include "FeederTask.h"

#include "GeantConfig.h"
#include "GeantRunManager.h"
#include "GeantScheduler.h"

#ifdef GEANT_TBB
#include "tbb/task_scheduler_init.h"
#endif

using namespace Geant;

FeederTask::FeederTask (Geant::GeantTaskData *td, bool starting): fTd(td), fStarting(starting) { }

FeederTask::~FeederTask () { }

tbb::task* FeederTask::execute ()
{
  GeantPropagator *prop = fTd->fPropagator;
  GeantRunManager *runmgr = prop->fRunMgr;
  int ninjected = 0;
  if (fStarting) {
//    while (runmgr->GetFedPropagator() != prop) {
//      int nb0 = runmgr->Feeder(fTd);
//      if (nb0 > 0) ninjected += nb0;
//    }
    if (!ninjected) prop->fWMgr->GetScheduler()->GarbageCollect(fTd, true);
  } else {    
//    ninjected = runmgr->Feeder(fTd);
  }
  
  tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
  TransportTask & transportTask = *new(cont.allocate_child()) TransportTask( fTd, fStarting );
  return & transportTask;
}
