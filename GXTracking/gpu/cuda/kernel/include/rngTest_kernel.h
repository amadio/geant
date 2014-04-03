#ifndef RNGTEST_KERNEL_H
#define RNGTEST_KERNEL_H

#include "GPTypeDef.h"
#include <curand_kernel.h>

extern void XORWOW_gpu(curandState *devStates,
		       int *result,
		       size_t nTrackSize,
		       int nBlocks,
		       int nThreads); 

extern void Mtgp32_gpu(curandStateMtgp32 *devStates,
		       int *result,
		       size_t nTrackSize,
		       int nBlocks,
		       int nThreads); 

extern void MRG32k3a_gpu(curandStateMRG32k3a *devStates,
			 int *result,
			 size_t nTrackSize,
			 int nBlocks,
			 int nThreads); 

extern void rand_cpu(int *result,
		     size_t nTrackSize);

#endif
