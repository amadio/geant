#ifndef TRACKING_GPU_H
#define TRACKING_GPU_H

#include "GPGeomManager.h"
#include "GPVPhysicalVolume.h"
#include "GXFieldMap.h"
#include "GXTrack.h"

#include "GPPhysicsTable.h"
#include <cmath>

#include <curand.h>
#include <curand_kernel.h>

extern void tracking_gpu(curandState *devStates, 
                         GPGeomManager *geomManager,
			 GXFieldMap *magMap,
			 GXTrack *track,  
			 GPPhysicsTable* eBrem_table, 
			 GPPhysicsTable* eIoni_table, 
			 GPPhysicsTable* msc_table, 
			 size_t nTrackSize,
			 int NBLOCKS,
			 int NTHREADS,
			 cudaStream_t stream);

extern void tracking_cpu( GPGeomManager *geomManager,
			 GXFieldMap *magMap, 
			 GXTrack *track, 
			 GPPhysicsTable* eBrem_table, 
			 GPPhysicsTable* eIoni_table, 
			 GPPhysicsTable* msc_table, 
			 size_t nTrackSize); 

extern void tracking_electron_gpu(curandState *devStates, 
                                  GPGeomManager *geomManager,
				  GXFieldMap *magMap,
				  GXTrack *track,  
				  GPPhysicsTable* eBrem_table, 
				  GPPhysicsTable* eIoni_table, 
				  GPPhysicsTable* msc_table, 
				  size_t nTrackSize,
				  int NBLOCKS,
				  int NTHREADS,
				  cudaStream_t stream);

extern void tracking_electron_cpu(GPGeomManager *geomManager,
				  GXFieldMap *magMap, 
				  GXTrack *track, 
				  GPPhysicsTable* eBrem_table, 
				  GPPhysicsTable* eIoni_table, 
				  GPPhysicsTable* msc_table, 
				  size_t nTrackSize); 

extern void tracking_photon_gpu(curandState *devStates,
				GPGeomManager *geomManager,
				GXFieldMap *magMap,
				GXTrack *track,  
				GPPhysicsTable* eBrem_table, 
				GPPhysicsTable* eIoni_table, 
				GPPhysicsTable* msc_table, 
				size_t nTrackSize,
				int NBLOCKS,
				int NTHREADS,
				cudaStream_t stream);

extern void tracking_photon_cpu(GPGeomManager *geomManager,
				GXFieldMap *magMap, 
				GXTrack *track, 
				GPPhysicsTable* eBrem_table, 
				GPPhysicsTable* eIoni_table, 
				GPPhysicsTable* msc_table, 
				size_t nTrackSize); 

extern void random_testing_gpu(curandState* devStates, double *result, 
                               int  blocksPerGrid, int threadsPerBlock,
                               cudaStream_t stream);

#endif
