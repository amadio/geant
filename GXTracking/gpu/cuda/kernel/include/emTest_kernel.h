#ifndef EMTEST_KERNEL_H
#define EMTEST_KERNEL_H

#include "GXTrack.h"
#include "GPPhysicsTable.h"
#include "GPGeomManager.h"
#include "GXFieldMap.h"

#include <curand.h>
#include <curand_kernel.h>

//the maximum number of secondaries per step
const G4int maxSecondaryPerStep = 2;

// electron_kernel wrapper
void electron_gpu(curandState *devStates,
		  GXTrack *track, size_t nTrackSize,
                  GPGeomManager *geomManager,
                  GXFieldMap *magMap,
		  GPPhysicsTable* physicsTable, 
		  GPPhysics2DVector* sbData,
		  GXTrack *secTrack,
		  G4int *stackSize,
		  G4int *offset,
		  G4int isStrack,
		  G4int runType,
		  int nBlocks,
		  int nThreads); 

void electron_cpu(GXTrack *trackIn, size_t nTrackSize,
                  GPGeomManager *geomManager,
                  GXFieldMap *magMap,
		  GPPhysicsTable* physicsTable,
		  GPPhysics2DVector* sbData,
		  GXTrack *secTrack,
		  G4int *stackSize,
		  G4int isStrack,
		  G4int runType);

// photon_kernel wrapper
void photon_gpu(curandState *devStates,
		GXTrack *track, size_t nTrackSize,
                GPGeomManager *geomManager,
                GXFieldMap *magMap,
		GPPhysicsTable* physicsTable,
		GXTrack *secTrack,
		G4int *stackSize,
		G4int *offset,
		G4int isStrack,
		G4int runType,
		int nBlocks,
		int nThreads); 

void photon_cpu(GXTrack *trackIn, size_t nTrackSize,
                GPGeomManager *geomManager,
                GXFieldMap *magMap,
		GPPhysicsTable* physicsTable,
		GXTrack *secTrack,
		G4int *stackSize,
		G4int isStrack,
		G4int runType);
#endif
