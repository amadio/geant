#ifndef RANDOM_KERNEL_H
#define RANDOM_KERNEL_H

#include <curand.h>
#include <curand_kernel.h>

extern void curand_setup_gpu(curandState* devStates, 
			     unsigned long seed,  
			     int NBLOCKS, 
			     int NTHREADS); 

extern void curand_setup_gpu(curandStateMtgp32* devStates, 
			     unsigned long seed,  
			     int NBLOCKS, 
			     int NTHREADS); 

extern void curand_setup_gpu(curandStateMRG32k3a* devStates, 
			     unsigned long seed,  
			     int NBLOCKS, 
			     int NTHREADS); 

#endif
