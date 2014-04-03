#ifndef GXMICManager_H
#define GXMICManager_H

#include "GXVCoprocessorManager.h"

class GXMICManager : public GXVCoprocessorManager
{
public:

  GXMICManager();
  ~GXMICManager();

  void Initialize();

  void AllocateDeviceMemory(/* GXTaskData* taskData_h */
			    GXVGeometry *geom,
			    GXFieldMap** fieldMap2D,
			    GXPhysicsTable* physicsTable,
			    GXPhysics2DVector* sbData);

  void DeallocateDeviceMemory();

  void UploadTaskData();

private:
  int numBlocks;
  int numThreads;

  GXVGeometry *geom_d;
  GXFieldMap* fieldMap_d;
  GXPhysicsTable* physicsTable_d;
  GXPhysics2DVector* sbData_d;

};

#endif
