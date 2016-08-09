#ifndef IOTASK
#define IOTASK

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

#endif //IOTASK
