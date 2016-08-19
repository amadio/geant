#ifndef GEANT_TBB_TRANSPORT_TASK
#define GEANT_TBB_TRANSPORT_TASK

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

#endif // GEANT_TBB_TRANSPORT_TASK
