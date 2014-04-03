#ifndef STEPPER_KERNEL_H
#define STEPPER_KERNEL_H

#include "GXFieldMap.h"
#include "GXTrack.h"

// rk4_kernel wrapper
extern void rk4_gpu(GXFieldMap *magMap,
		    GXTrack *track, 
		    size_t nTrackSize,
		    int nBlocks,
		    int nThreads); 

extern void rk4_cpu(GXFieldMap *magMap,
		    GXTrack *trackIn,
		    size_t nTrackSize);

// rkf45_kernel wrapper
extern void rkf45_gpu(GXFieldMap *magMap,
		      GXTrack *track, 
		      size_t nTrackSize,
		      int nBlocks,
		      int nThreads); 

extern void rkf45_cpu(GXFieldMap *magMap,
		      GXTrack *trackIn,
		      size_t nTrackSize);

// nrk4_kernel wrapper
extern void nrk4_gpu(GXFieldMap *magMap,
		     GXTrack *track, 
		     size_t nTrackSize,
		     int nBlocks,
		     int nThreads); 

extern void nrk4_cpu(GXFieldMap *magMap,
		     GXTrack *trackIn,
		     size_t nTrackSize);

#endif
