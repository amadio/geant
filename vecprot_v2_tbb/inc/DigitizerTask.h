#ifndef DIGITIZERTASK
#define DIGITIZERTASK

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

#endif //DIGITIZERTASK
