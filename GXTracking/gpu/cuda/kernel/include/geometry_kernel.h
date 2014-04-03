#ifndef GEOMETRY_GPU_H
#define GEOMETRY_GPU_H

// geometry_kernel wrapper

#include "GPGeomManager.h"
#include "GPVPhysicalVolume.h"
#include "GPNavigator.h"
#include "GXFieldMap.h"
#include "GXTrack.h"

extern void navigator_gpu(GPGeomManager *geomManager,
                          GXFieldMap *magMap,
                          GXTrack *track,  
                          size_t nTrackSize,
                          int NBLOCKS,
                          int NTHREADS);

extern void navigator_cpu(GPGeomManager *geomManager,
                          GXFieldMap *magMap, 
                          GXTrack *track, 
                          size_t nTrackSize); 

extern void mllocator_gpu(GPGeomManager *geomManager,
			  GXFieldMap *magMap,
			  GXTrack *track,  
			  size_t nTrackSize,
			  int NBLOCKS,
			  int NTHREADS);

extern void mllocator_cpu(GPGeomManager *geomManager,
			  GXFieldMap *magMap, 
			  GXTrack *track, 
			  size_t nTrackSize); 


#endif
