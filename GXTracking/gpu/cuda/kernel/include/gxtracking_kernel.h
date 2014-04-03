#ifndef GXTRACKING_KERNEL_H
#define GXTRACKING_KERNEL_H

#include "GXTrack.h"
#include "GXPhysicsTable.h"
#include "GXPhysics2DVector.h"
#include "GPGeomManager.h"
#include "GXFieldMap.h"

#include <curand.h>
#include <curand_kernel.h>

//the maximum number of secondaries per step
//const G4int maxSecondaryPerStep = 2;

// elec_GPIL_kernel wrapper
void elec_GPIL_gpu(curandState *devStates,
		   GXTrack *track, 
		   GXTrackLiason *liason, 
		   size_t nTrackSize,
		   GPGeomManager *geomManager,
		   GXFieldMap *magMap,
		   GXPhysicsTable* physicsTable, 
		   GXPhysics2DVector* sbData,
		   GXTrack *secTrack,
		   G4int *stackSize,
		   G4int numStep,
		   int nBlocks,
		   int nThreads,
		   cudaStream_t stream); 

void elec_GPIL_cpu(GXTrack *trackIn, 
		   GXTrackLiason *liason, 
		   size_t nTrackSize,
		   GPGeomManager *geomManager,
		   GXFieldMap *magMap,
		   GXPhysicsTable* physicsTable,
		   GXPhysics2DVector* sbData,
		   GXTrack *secTrack,
		   G4int *stackSize,
		   G4int numStep);

#endif
