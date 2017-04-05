#include "TaskMgrTBB.h"

#include "tbb/task_scheduler_init.h"
#include "tbb/task.h"
#include "InitialTask.h"
#include "ThreadData.h"

using namespace Geant;

//______________________________________________________________________________
bool TaskMgrTBB::Initialize(int nthreads, GeantPropagator *prop)
{
// Initialize the TBB task system

  // Start Task scheduler and tasks
  tbb::task_list tlist;
  // spawn transport tasks
  for (int i = 0; i < nthreads; i++)
    tlist.push_back(*new (tbb::task::allocate_root()) InitialTask(prop));

  tbb::task::spawn_root_and_wait(tlist);
  return true;
}

//______________________________________________________________________________
void TaskMgrTBB::Finalize()
{
// Final actions for the TBB task system
  // Cleanup task data
  delete ThreadData::Instance(0);
}
