#ifndef TRANSPORTTASK
#define TRANSPORTTASK

#include "GeantTaskData.h"
#include "tbb/task.h"

class TransportTask : public tbb::task
{
private:
  Geant::GeantTaskData *fTd;
  int fNbaskets;

public:
  TransportTask (Geant::GeantTaskData *td, int nbaskets);
  ~TransportTask ();

  tbb::task* execute ();

};

#endif // TRANSPORTTASK
