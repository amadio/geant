#include "GUCurand.h"

#include <iostream>

namespace vecphys {
inline namespace cuda {

__global__
void GUCurand_Init_Kernel(Random_t *randomStates, unsigned long seed) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &randomStates[tid]);
}

} // end of cuda namespace

bool GUCurand_Init(Random_t *randomStates, unsigned long seed, 
                   int blocksPerGrid, int threadsPerBlock) 
{
  int kstatus = 0;

  vecphys::cuda::GUCurand_Init_Kernel<<<blocksPerGrid,threadsPerBlock>>>
                                                     (randomStates,seed);

  cudaError err = cudaGetLastError();
  if ( cudaSuccess != err ) {
    fprintf(stderr,"GUCurand_Init cudaCheckError() failed at %s : %i : %s\n",
	    __FILE__, __LINE__, cudaGetErrorString(err));
    exit(-1);
  }
  
  kstatus = cudaThreadSynchronize();
  if (kstatus) std::cout << "GUCurand_Init status = " << kstatus << "\n";
  
  return true;
}

} // end of vecphys namespace
