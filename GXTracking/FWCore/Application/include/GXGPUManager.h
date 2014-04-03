#ifndef GXGPUManager_H
#define GXGPUManager_H

#include "GXVCoprocessorManager.h"
#include "GXTrackHandler.h"

class GXGPUManager : public GXVCoprocessorManager
{
public:

  GXGPUManager(int NBlocks = 1, int NThreads = 1);
  ~GXGPUManager();

  void Initialize();
  void AllocateDeviceMemory(/* GXTaskData* taskData_h */
			    GXVGeometry *geom,
			    GXFieldMap** fieldMap2D,
			    GXPhysicsTable* physicsTable,
			    GXPhysics2DVector* sbData);
  void DeallocateDeviceMemory();
  void UploadTaskData();
  void DownloadTaskData();
  void DeallocateTaskData();
  void LaunchTask();

  void SetBlockThread(int nblocks, int nthreads);
  void SetNumberOfSteps(int nsteps);
  void SetLimitHeapSize(size_t bytes);

  void StartTimer();
  float StopTimer();

  inline GXTrackHandler* GetTrackHandler() { return fTrackHandler; }; 

private:

  int fNTracks;
  int fNSteps;
  curandState* fRandomStates;
  GXTrackHandler* fTrackHandler;

  int fNBlocks;
  int fNThreads;
  size_t fLimitHeapSize;

  GXVGeometry::byte *geom_d;
  GXFieldMap* fieldMap_d;
  GXPhysicsTable* physicsTable_d;
  GXPhysics2DVector* sbData_d;
  GXTrack* tracks_d;

  cudaStream_t stream;
  cudaEvent_t start;
  cudaEvent_t stop;

  //fixed size secondary stack and counters
  GXTrack *secTracks_d;  
  G4int *stackSize_d;

//--->@@@G4FWP - these are temporary members for CPU tests
public:
  void LaunchCPUTask();
  void PrintPerformance(int taskId);
  void SetPerformanceFlag(bool flag) { fPerformance = flag; };

  bool fPerformance;
  float fElapsedTimeH2D;
  float fElapsedTimeD2H;
  float fElapsedTimeGPU;
  float fElapsedTimeCPU;

private:
  GXVGeometry::byte *geom_h;
  GXFieldMap* fieldMap_h;
  GXPhysicsTable* physicsTable_h;
  GXPhysics2DVector* sbData_h;
  GXTrack* tracks_h;
//<---@@@G4FWP

};

#endif
