#include "TaskMgrTBB.h"

#include "tbb/task_scheduler_init.h"
#include "tbb/task.h"
#include "InitialTask.h"

//______________________________________________________________________________
bool TaskMgrTBB::Initialize(int nthreads)
{
// Initialize the TBB task system

  // Start Task scheduler and tasks
  tbb::task_list tlist;
  tbb::task &cont = *new (tbb::task::allocate_root()) tbb::empty_task();
  // spawn transport tasks
  for (int i = 0; i < nthreads; i++)
    tlist.push_back(*new (tbb::task::allocate_root()) InitialTask());

  tbb::task::spawn_root_and_wait(tlist);
  return true;
}
