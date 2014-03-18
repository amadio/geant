#include <iostream>

#include "GPTypeDef.h"
#include <curand_kernel.h>
#include "GPRandom.h"

//-----------------------------------------------------------------------------
//  XORWOW_kernel
//-----------------------------------------------------------------------------
GLOBALFUNC 
void XORWOW_kernel(curandState *devStates, int *result, size_t nTracks) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  int count = 0;
  G4double x;

  curandState localState = devStates[tid];

  for(int i=0 ; i < nTracks ; i++) {
    x = curand_uniform(&localState);
    if ( x < 0.5 ) count++;  
  }
  devStates[tid] = localState;

  result[tid] += count;
}

//-----------------------------------------------------------------------------
//  cuda wrapper for XORWOW_kernel
//-----------------------------------------------------------------------------
void XORWOW_gpu(curandState *devStates, int *result, size_t nTrackSize,
		int NBLOCKS, int NTHREADS) 
{
  int kstatus = 0;

  XORWOW_kernel<<<NBLOCKS,NTHREADS>>>(devStates,result,nTrackSize);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "XORWOW_gpu status = " << kstatus << "\n";
}

//-----------------------------------------------------------------------------
//  curandStateMtgp32_kernel
//-----------------------------------------------------------------------------
GLOBALFUNC 
void Mtgp32_kernel(curandStateMtgp32 *devStates, int *result, size_t nTracks) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  int count = 0;
  G4double x;

  for(int i=0 ; i < nTracks ; i++) {
    x = curand_uniform(&devStates[blockIdx.x]);
    if ( x <0.5 ) count++;  
  }
  result[tid] += count;
}

//-----------------------------------------------------------------------------
//  cuda wrapper for Mtgp32_kernel
//-----------------------------------------------------------------------------
void Mtgp32_gpu(curandStateMtgp32 *devStates, int *result, size_t nTrackSize,
		int NBLOCKS, int NTHREADS) 
{
  int kstatus = 0;

  Mtgp32_kernel<<<NBLOCKS,NTHREADS>>>(devStates,result,nTrackSize);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "Mtgp32_gpu status = " << kstatus << "\n";
}

//-----------------------------------------------------------------------------
//  MRG32k3a_kernel
//-----------------------------------------------------------------------------
GLOBALFUNC 
void MRG32k3a_kernel(curandStateMRG32k3a *devStates, int *result, 
		     size_t nTracks) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  int count = 0;
  G4double x;

  curandStateMRG32k3a localState = devStates[tid];

  for(int i=0 ; i < nTracks ; i++) {
    x = curand_uniform(&localState);
    if ( x < 0.5 ) count++;  
  }
  devStates[tid] = localState;

  result[tid] += count;
}

//-----------------------------------------------------------------------------
//  cuda wrapper for MRG32k3a_kernel
//-----------------------------------------------------------------------------
void MRG32k3a_gpu(curandStateMRG32k3a *devStates, int *result, 
		  size_t nTrackSize, int NBLOCKS, int NTHREADS) 
{
  int kstatus = 0;

  MRG32k3a_kernel<<<NBLOCKS,NTHREADS>>>(devStates,result,nTrackSize);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "MRG32k3a_gpu status = " << kstatus << "\n";
}

void rand_cpu(int *result, size_t nTracks) 
{
  curandState *devStates;
  int count = 0;
  G4double x;

  for(int i=0 ; i < nTracks ; i++) {
    x = rand_wrapper(devStates,0);
    if ( x < 0.5 ) count++;  
  }
  result[0] = count;
}

#include "random_kernel.cu"
