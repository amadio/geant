#ifndef GXVCoprocessorManager_H
#define GXVCoprocessorManager_H

#include "GXVGeometry.h"
#include "GXFieldMap.h"
#include "GXPhysicsTable.h"
#include "GXPhysics2DVector.h"

struct GXTaskData {
  GXVGeometry* geom;
  GXFieldMap** fieldMap;
  GXPhysicsTable* physicsTable;
  GXPhysics2DVector* sbData;
};

class GXVCoprocessorManager
{
public:
	
  GXVCoprocessorManager(char* name) { coprocessorName = name; }
  virtual ~GXVCoprocessorManager() {}
  virtual char* GetName() { return coprocessorName ; };

  virtual void Initialize() = 0;  
  virtual void AllocateDeviceMemory(/* GXTaskData* taskData */
				    GXVGeometry *geom,
				    GXFieldMap** fieldMap,
				    GXPhysicsTable* physicsTable,
				    GXPhysics2DVector* sbData) = 0;

  virtual void DeallocateDeviceMemory() = 0;


  virtual void UploadTaskData() = 0;

private:
  char* coprocessorName;
};

#endif
