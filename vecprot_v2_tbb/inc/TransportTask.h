#ifndef TRANSPORTTASK
#define TRANSPORTTASK

#include "GeantTaskData.h"
#include "tbb/task.h"

class TransportTask : public tbb::task
{
private:
  Geant::GeantTaskData *fTd;

public:
  TransportTask (Geant::GeantTaskData *td);
  ~TransportTask ();

  tbb::task* execute ();

};

#endif // TRANSPORTTASK
