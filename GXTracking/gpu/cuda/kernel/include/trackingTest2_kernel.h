#ifndef TRACKINGTEST2_KERNEL_H
#define TRACKINGTEST2_KERNEL_H

#include "GXTrack.h"
#include "GPPhysicsTable.h"
#include "GPGeomManager.h"
#include "GXFieldMap.h"

#include <curand.h>
#include <curand_kernel.h>

class GPPhysics2DVector;
class GXTrackLiason;

//the maximum number of secondaries per step
const G4int maxSecondaryPerStep = 2;

// elec_kernel wrapper
void elec_gpu(curandState *devStates,
              GXTrack *track,
              size_t nTrackSize,
              GPGeomManager *geomManager,
              GXFieldMap *magMap,
              GPPhysicsTable* physicsTable,
              GPPhysics2DVector* sbData,
              GXTrack *secTrack,
              G4int *stackSize,
              G4int isStrack,
              G4int runType,
              int nBlocks,
              int nThreads,
              cudaStream_t stream);

void elec_gpu(curandState *devStates,
              GXTrack *track,
              int *logVolumeIndices,
              int *physVolumeIndices,
              size_t nTrackSize,
              GPGeomManager *geomManager,
              GXFieldMap *magMap,
              GPPhysicsTable* physicsTable,
              GPPhysics2DVector* sbData,
              GXTrack *secTrack,
              G4int *stackSize,
              G4int isStrack,
              G4int runType,
              int nBlocks,
              int nThreads,
              cudaStream_t stream);

void elec_cpu(GXTrack *trackIn, 
              size_t nTrackSize,
              GPGeomManager *geomManager,
              GXFieldMap *magMap,
              GPPhysicsTable* physicsTable,
              GPPhysics2DVector* sbData,
              GXTrack *secTrack,
              G4int *stackSize,
              G4int isStrack,
              G4int runType);

// elec_GPIL_kernel wrapper
void elec_GPIL_gpu(curandState *devStates,
                   GXTrack *track,
                   GXTrackLiason *liason,
                   size_t nTrackSize,
                   GPGeomManager *geomManager,
                   GXFieldMap *magMap,
                   GPPhysicsTable* physicsTable,
                   GPPhysics2DVector* sbData,
                   GXTrack *secTrack,
                   G4int *stackSize,
                   G4int isStrack,
                   G4int runType,
                   int nBlocks,
                   int nThreads,
                   cudaStream_t stream);

void elec_GPIL_gpu(curandState *devStates,
                   GXTrack *track,
                   GXTrackLiason *liason,
                   int *logVolumeIndices,
                   int *physVolumeIndices,
                   size_t nTrackSize,
                   GPGeomManager *geomManager,
                   GXFieldMap *magMap,
                   GPPhysicsTable* physicsTable,
                   GPPhysics2DVector* sbData,
                   GXTrack *secTrack,
                   G4int *stackSize,
                   G4int isStrack,
                   G4int runType,
                   int nBlocks,
                   int nThreads,
                   cudaStream_t stream);

void elec_GPIL_cpu(GXTrack *trackIn, 
                   GXTrackLiason *liason,
                   size_t nTrackSize,
                   GPGeomManager *geomManager,
                   GXFieldMap *magMap,
                   GPPhysicsTable* physicsTable,
                   GPPhysics2DVector* sbData,
                   GXTrack *secTrack,
                   G4int *stackSize,
                   G4int isStrack,
                   G4int runType);

// elec_doit_kernel wrapper
void elec_doit_gpu(curandState *devStates,
                   GXTrack *track,
                   GXTrackLiason *liason,
                   size_t nTrackSize,
                   GPGeomManager *geomManager,
                   GXFieldMap *magMap,
                   GPPhysicsTable* physicsTable,
                   GPPhysics2DVector* sbData,
                   GXTrack *secTrack,
                   G4int *stackSize,
                   G4int isStrack,
                   G4int runType,
                   int nBlocks,
                   int nThreads,
                   cudaStream_t stream);

void elec_doit_cpu(GXTrack *trackIn, 
                   GXTrackLiason *liason,
                   size_t nTrackSize,
                   GPGeomManager *geomManager,
                   GXFieldMap *magMap,
                   GPPhysicsTable* physicsTable,
                   GPPhysics2DVector* sbData,
                   GXTrack *secTrack,
                   G4int *stackSize,
                   G4int isStrack,
                   G4int runType);

#endif
