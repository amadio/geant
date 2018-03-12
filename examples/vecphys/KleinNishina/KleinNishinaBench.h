#ifndef GEANTV_KLEINNISHINABENCH_H
#define GEANTV_KLEINNISHINABENCH_H

#include "KleinNishinaTestCommon.h"

void PreparePrimaries(std::vector<LightTrack> &output);
TaskData *PrepareTaskData();
void CleanTaskData(TaskData *Td);

#endif // GEANTV_KLEINNISHINABENCH_H
