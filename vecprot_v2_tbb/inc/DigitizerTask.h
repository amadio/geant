#ifndef GEANT_TBB_DIGITIZER_TASK
#define GEANT_TBB_DIGITIZER_TASK

#include "GeantTaskData.h"
#include "tbb/task.h"

class DigitizerTask : public tbb::task
{
private:
  Geant::GeantTaskData *fTd;

public:
  DigitizerTask (Geant::GeantTaskData *td);
  ~DigitizerTask ();

  tbb::task* execute ();

};

#endif //GEANT_TBB_DIGITIZER_TASK
