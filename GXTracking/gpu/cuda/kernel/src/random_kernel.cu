#include <cuda.h>
#include <curand_kernel.h>

#include <iostream>
#include "GPTypeDef.h"

#ifdef __CUDACC__
GLOBALFUNC void curand_setup_kernel(curandState *devStates, unsigned long seed) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &devStates[tid]);
}

GLOBALFUNC void curand_setup_kernel(curandStateMRG32k3a *devStates, unsigned long seed) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &devStates[tid]);

}

void curand_setup_gpu(curandState *devStates, unsigned long seed,  int NBLOCKS, int NTHREADS) {

  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid = NBLOCKS;

  curand_setup_kernel<<< blocksPerGrid, threadsPerBlock >>> (devStates,seed);

  kstatus = cudaThreadSynchronize();
  if (kstatus)
    std::cout << "cuda_setup_kernel status = " << kstatus << "\n";

}

void curand_setup_gpu(curandStateMRG32k3a *devStates, unsigned long seed,  int NBLOCKS, int NTHREADS) {

  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid = NBLOCKS;

  curand_setup_kernel<<< blocksPerGrid, threadsPerBlock >>> (devStates,seed);

  kstatus = cudaThreadSynchronize();
  if (kstatus)
    std::cout << "cuda_setup_kernel status = " << kstatus << "\n";

}
#endif

