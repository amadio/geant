#ifndef GEANTV_COMMON_H
#define GEANTV_COMMON_H

#include "Geant/PhysicsData.h"
using geant::TaskData;

TaskData *PrepareTaskData()
{
  auto PhysData    = new geantphysics::PhysicsData();
  auto Td          = new TaskData(1, 2048);
  Td->fPhysicsData = PhysData;
  return Td;
}

void CleanTaskData(TaskData *Td)
{
  delete Td->fPhysicsData;
  delete Td;
}

#endif // GEANTV_COMMON_H
