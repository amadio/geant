#ifndef TRACKINGTEST_KERNEL_H
#define TRACKINGTEST_KERNEL_H

#include "GXTrack.h"
#include "GPPhysicsTable.h"
#include "GPGeomManager.h"
#include "GXFieldMap.h"

#include <curand.h>
#include <curand_kernel.h>

//the maximum number of secondaries per step
const G4int maxSecondaryPerStep = 2;

// brem_kernel wrapper
void brem_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
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

void brem_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType);

// ioni_kernel wrapper
void ioni_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
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

void ioni_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable,
	      GPPhysics2DVector* sbData,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType);

// msc_kernel wrapper
void msc_gpu(curandState *devStates,
	     GXTrack *track, size_t nTrackSize,
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

void msc_cpu(GXTrack *trackIn, size_t nTrackSize,
	     GPGeomManager *geomManager,
	     GXFieldMap *magMap,
	     GPPhysicsTable* physicsTable,
	     GPPhysics2DVector* sbData,
	     GXTrack *secTrack,
	     G4int *stackSize,
	     G4int isStrack,
	     G4int runType);

// trans_kernel wrapper
void trans_gpu(curandState *devStates,
	       GXTrack *track, size_t nTrackSize,
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

void trans_cpu(GXTrack *trackIn, size_t nTrackSize,
	       GPGeomManager *geomManager,
	       GXFieldMap *magMap,
	       GPPhysicsTable* physicsTable,
	       GPPhysics2DVector* sbData,
	       GXTrack *secTrack,
	       G4int *stackSize,
	       G4int isStrack,
	       G4int runType);

// compt_kernel wrapper
void compt_gpu(curandState *devStates,
	       GXTrack *track, size_t nTrackSize,
	       GPGeomManager *geomManager,
	       GXFieldMap *magMap,
	       GPPhysicsTable* physicsTable,
	       GXTrack *secTrack,
	       G4int *stackSize,
	       G4int isStrack,
	       G4int runType,
	       int nBlocks,
	       int nThreads,
	       cudaStream_t stream);

void compt_cpu(GXTrack *trackIn, size_t nTrackSize,
	       GPGeomManager *geomManager,
	       GXFieldMap *magMap,
	       GPPhysicsTable* physicsTable,
	       GXTrack *secTrack,
	       G4int *stackSize,
	       G4int isStrack,
	       G4int runType);

// conv_kernel wrapper
void conv_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType,
	      int nBlocks,
	      int nThreads,
	      cudaStream_t stream);

void conv_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType);

// phot_kernel wrapper
void phot_gpu(curandState *devStates,
	      GXTrack *track, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType,
	      int nBlocks,
	      int nThreads,
	      cudaStream_t stream);

void phot_cpu(GXTrack *trackIn, size_t nTrackSize,
	      GPGeomManager *geomManager,
	      GXFieldMap *magMap,
	      GPPhysicsTable* physicsTable,
	      GXTrack *secTrack,
	      G4int *stackSize,
	      G4int isStrack,
	      G4int runType);
#endif
