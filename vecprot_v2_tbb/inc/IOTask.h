#ifndef GEANT_TBB_IO_TASK
#define GEANT_TBB_IO_TASK

#include "GeantTaskData.h"
#include "tbb/task.h"

class IOTask : public tbb::task
{
private:
  Geant::GeantTaskData *fTd;

public:
  IOTask (Geant::GeantTaskData *td);
  ~IOTask ();

  tbb::task* execute ();

};

#endif //GEANT_TBB_IO_TASK
