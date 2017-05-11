#ifndef RNG_BENCHMARKER_GPU_H
#define RNG_BENCHMARKER_GPU_H 1

#include "rng/MRG32k3a.h"
#include "rng/Threefry.h"
#include "rng/Philox.h"

#include <curand_kernel.h>

namespace vecphys {

  // Cuda Wrapper
  void CudaMRG32k3a(MRG32k3a_t<ScalarBackend> *devStates, double *result, 
		    int nsample, int blocksPerGrid, int threadsPerBlock);

  void CudaThreefry(Threefry_t<ScalarBackend> *devStates, double *result, 
		    int nsample, int blocksPerGrid, int threadsPerBlock);

  void CudaPhilox(Philox_t<ScalarBackend> *devStates, double *result, 
		  int nsample, int blocksPerGrid, int threadsPerBlock);

  // Curand Test 
  void CurandMRG32k3a(curandStateMRG32k3a *devStates, double *result,
		      int nsample, int blocksPerGrid, int threadsPerBlock);

  void curand_setup_gpu(curandStateMRG32k3a* devStates, unsigned long seed,  
                        int blocksPerGrid, int threadsPerBlock);

  void CurandPhilox(curandStatePhilox4_32_10_t *devStates, double *result,
                    int nsample, int blocksPerGrid, int threadsPerBlock);

  void curand_setup_gpu(curandStatePhilox4_32_10_t* devStates, unsigned long seed,  
                        int blocksPerGrid, int threadsPerBlock);

} // end namespace vecphys

#endif
