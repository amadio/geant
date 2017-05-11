#include "rng/MRG32k3a.h"
#include "rng/Threefry.h"
#include "rng/Philox.h"

namespace vecphys {
inline namespace cuda {

__global__
void KernelMRG32k3a(MRG32k3a_t<ScalarBackend>* devStates, double *result, int nsample) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = tid;

  MRG32k3a<ScalarBackend> rng(0);

  __shared__ double sum[26*192];
  double tmp = 0;

  while (tid < nsample) {
    tmp += rng.Uniform<ScalarBackend>(&devStates[sid]);
    tid += blockDim.x * gridDim.x;
  }

  sum[sid] = tmp;

  __syncthreads();

  //do reduction on CPU
  result[sid] = sum[sid];
}

__global__
void KernelThreefry(Threefry_t<ScalarBackend>* devStates, double *result, int nsample) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = tid;

  Threefry<ScalarBackend> rng(0);
 
  __shared__ double sum[26*192];
  double tmp = 0;

  while (tid < nsample) {
    tmp += rng.Uniform<ScalarBackend>(&devStates[sid]);
    tid += blockDim.x * gridDim.x;
  }

  sum[sid] = tmp;

  __syncthreads();

  //do reduction on CPU
  result[sid] = sum[sid];
}

__global__
void KernelPhilox(Philox_t<ScalarBackend>* devStates, double *result, int nsample) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = tid;

  Philox<ScalarBackend> rng(0);
 
  __shared__ double sum[26*192];
  double tmp = 0;

  while (tid < nsample) {
    tmp += rng.Uniform<ScalarBackend>(&devStates[sid]);
    tid += blockDim.x * gridDim.x;
  }

  sum[sid] = tmp;

  __syncthreads();

  //do reduction on CPU
  result[sid] = sum[sid];
}

//-----------------------------------------------------------------------------
//  Curand MRG32k3a
//-----------------------------------------------------------------------------
__global__
void KernelCurandMRG32k3a(curandStateMRG32k3a *devStates, double *result, int nsample)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = tid;

  curandStateMRG32k3a localState = devStates[sid];

  __shared__ double sum[26*192];
  double tmp = 0;

  while (tid < nsample) {
    tmp += curand_uniform_double(&localState);
    //    tmp += devStates[sid].s1[0];
    tid += blockDim.x * gridDim.x;
  }

  devStates[sid] = localState;

  sum[sid] = tmp;

  __syncthreads();

  //do reduction on CPU
  result[sid] = sum[sid];
}

__global__
void curand_setup_kernel(curandStateMRG32k3a *devStates, unsigned long seed) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &devStates[tid]);
}

//-----------------------------------------------------------------------------
//  Curand Philox
//-----------------------------------------------------------------------------
__global__
void KernelCurandPhilox(curandStatePhilox4_32_10_t *devStates, double *result, int nsample)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = tid;

  curandStatePhilox4_32_10_t localState = devStates[tid];

  __shared__ double sum[26*192];
  double tmp = 0;

  while (tid < nsample) {
    tmp += curand_uniform_double(&localState);
    tid += blockDim.x * gridDim.x;
  }

  devStates[sid] = localState;

  sum[sid] = tmp;

  __syncthreads();

  //do reduction on CPU
  result[sid] = sum[sid];
}


__global__
void curand_setup_kernel(curandStatePhilox4_32_10_t *devStates, unsigned long seed) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &devStates[tid]);
}


} // end namespace cuda

// Cuda wrapper

void CudaMRG32k3a(MRG32k3a_t<ScalarBackend> *devStates,
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
  KernelMRG32k3a<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);
}

void CudaThreefry(Threefry_t<ScalarBackend> *devStates,
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
   KernelThreefry<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);
}

void CudaPhilox(Philox_t<ScalarBackend> *devStates,
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
   KernelPhilox<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for Curand MRG32k3a Kernel
//-----------------------------------------------------------------------------
void CurandMRG32k3a(curandStateMRG32k3a *devStates, 
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
  int kstatus = 0;

  KernelCurandMRG32k3a<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "MRG32k3a_gpu status = " << kstatus << "\n";
}


void curand_setup_gpu(curandStateMRG32k3a *devStates, unsigned long seed,  int NBLOCKS, int NTHREADS) {

  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid = NBLOCKS;

  curand_setup_kernel<<< blocksPerGrid, threadsPerBlock >>> (devStates,seed);

  kstatus = cudaThreadSynchronize();
  if (kstatus)
    std::cout << "MRG32k3a: cuda_setup_kernel status = " << kstatus << "\n";

}

//-----------------------------------------------------------------------------
//  cuda wrapper for Curand Philox4_32_10 Kernel
//-----------------------------------------------------------------------------
void CurandPhilox(curandStatePhilox4_32_10_t *devStates, 
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
  int kstatus = 0;

  KernelCurandPhilox<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "CurandPhilox status = " << kstatus << "\n";
}


void curand_setup_gpu(curandStatePhilox4_32_10_t *devStates, unsigned long seed,  int NBLOCKS, int NTHREADS) {

  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid = NBLOCKS;

  curand_setup_kernel<<< blocksPerGrid, threadsPerBlock >>> (devStates,seed);

  kstatus = cudaThreadSynchronize();
  if (kstatus)
    std::cout << "Philox: cuda_setup_kernel status = " << kstatus << "\n";

}

} // end namespace vecphys
