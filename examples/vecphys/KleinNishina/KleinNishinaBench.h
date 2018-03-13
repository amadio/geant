#ifndef GEANTV_KLEINNISHINABENCH_H
#define GEANTV_KLEINNISHINABENCH_H

#include "KleinNishinaTestCommon.h"

void PreparePrimaries(std::vector<LightTrack> &output, int N);
void PreparePrimaries(LightTrack_v &outputm, int N);
TaskData *PrepareTaskData();
void CleanTaskData(TaskData *Td);

#endif // GEANTV_KLEINNISHINABENCH_H
