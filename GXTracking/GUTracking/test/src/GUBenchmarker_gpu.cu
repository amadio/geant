#include "base/Stopwatch.h"

#include "GUComptonKleinNishina.h"
#include "GVComptonKleinNishina.h"
#include "GUConversionBetheHeitler.h"
#include "GUPhotoElectronSauterGavrila.h"
#include "GUMollerBhabha.h"
#include "GUSeltzerBerger.h"

#ifdef VECPHYS_ROOT
#include "GUHistogram.h"
#endif

namespace vecphys {

inline namespace cuda {

__global__
void KernelKleinNishina(Random_t* devStates,
                        GUAliasTableManager** table,
                        Physics2DVector* sbData,
			int nTrackSize, 
                        GUTrack* itrack, 
			int* targetElements, 
			GUTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GUAliasSampler sampler(devStates,tid,1.e-8,1.e+3,100,100,table[kKleinNishina]);
  GUComptonKleinNishina model(devStates,tid,&sampler);

  while (tid < nTrackSize) {
    model.Interact<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__
void KernelVKleinNishina(Random_t* devStates,
                         GUAliasTableManager** table,
                         Physics2DVector* sbData,
			 int nTrackSize, 
                         GUTrack* itrack, 
			 int* targetElements, 
			 GUTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GUAliasSampler sampler(devStates,tid,1.e-8,1.e+3,100,100,table[kVKleinNishina]);
  GVComptonKleinNishina model(devStates,tid,&sampler);

  while (tid < nTrackSize) {
    model.Interact<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__
void KernelBetheHeitler(Random_t* devStates,
                        GUAliasTableManager** table,
                        Physics2DVector* sbData,
			int nTrackSize, 
                        GUTrack* itrack, 
			int* targetElements, 
			GUTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GUAliasSampler sampler(devStates,tid,1.e-8,1.e+3,100,100,table[kBetheHeitler]);
  GUConversionBetheHeitler model(devStates,tid,&sampler);

  while (tid < nTrackSize) {
    model.Interact<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__
void KernelSauterGavrila(Random_t* devStates,
                         GUAliasTableManager** table,
                         Physics2DVector* sbData,
			 int nTrackSize, 
                         GUTrack* itrack, 
			 int* targetElements, 
			 GUTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GUAliasSampler sampler(devStates,tid,1.e-8,1.e+3,100,100,table[kSauterGavrila]);
  GUPhotoElectronSauterGavrila model(devStates,tid,&sampler);

  while (tid < nTrackSize) {
    model.Interact<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__
void KernelMollerBhabha(Random_t* devStates,
                        GUAliasTableManager** table,
                        Physics2DVector* sbData,
			int nTrackSize, 
                        GUTrack* itrack, 
			int* targetElements, 
			GUTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GUAliasSampler sampler(devStates,tid,1.e-8,1.e+3,100,100,table[kMollerBhabha]);
  GUMollerBhabha model(devStates,tid,&sampler);

  while (tid < nTrackSize) {
    model.Interact<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__
void KernelSeltzerBerger(Random_t* devStates,
                         GUAliasTableManager** table,
                         Physics2DVector* sbData,
			 int nTrackSize, 
                         GUTrack* itrack, 
			 int* targetElements, 
			 GUTrack* otrack)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GUAliasSampler sampler(devStates,tid,1.e-8,1.e+3,100,100,table[kSeltzerBerger]);
  GUSeltzerBerger model(devStates,tid,&sampler,sbData);

  while (tid < nTrackSize) {
    model.Interact<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

} // end namespace cuda

// Cuda 

Precision CudaKleinNishina(int blocksPerGrid, 
                           int threadsPerBlock,
                           Random_t* devStates,
                           GUAliasTableManager** table,
                           Physics2DVector* sbData,
			   int nTrackSize, 
                           GUTrack* itrack, 
			   int* targetElements, 
			   GUTrack* otrack)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelKleinNishina<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack);
  cudaDeviceSynchronize();

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision CudaVKleinNishina(int blocksPerGrid, 
                            int threadsPerBlock,
                            Random_t* devStates,
                            GUAliasTableManager** table,
                            Physics2DVector* sbData,
			    int nTrackSize, 
                            GUTrack* itrack, 
			    int* targetElements, 
			    GUTrack* otrack)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelVKleinNishina<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack);
  cudaDeviceSynchronize();

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision CudaBetheHeitler(int blocksPerGrid, 
                           int threadsPerBlock,
                           Random_t* devStates,
                           GUAliasTableManager** table,
                           Physics2DVector* sbData,
			   int nTrackSize, 
                           GUTrack* itrack, 
			   int* targetElements, 
			   GUTrack* otrack)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelBetheHeitler<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack);
  cudaDeviceSynchronize();

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision CudaSauterGavrila(int blocksPerGrid, 
                            int threadsPerBlock,
                            Random_t* devStates,
                            GUAliasTableManager** table,
                            Physics2DVector* sbData,
			    int nTrackSize, 
                            GUTrack* itrack, 
			    int* targetElements, 
			    GUTrack* otrack)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelSauterGavrila<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack);
  cudaDeviceSynchronize();

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision CudaMollerBhabha(int blocksPerGrid, 
                           int threadsPerBlock,
                           Random_t* devStates,
                           GUAliasTableManager** table,
                           Physics2DVector* sbData,
			   int nTrackSize, 
                           GUTrack* itrack, 
			   int* targetElements, 
			   GUTrack* otrack)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelMollerBhabha<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack);
  cudaDeviceSynchronize();

  elapsedTime = timer.Stop();

  return elapsedTime;
}

Precision CudaSeltzerBerger(int blocksPerGrid, 
                            int threadsPerBlock,
                            Random_t* devStates,
                            GUAliasTableManager** table,
                            Physics2DVector* sbData,
			    int nTrackSize, 
                            GUTrack* itrack, 
			    int* targetElements, 
			    GUTrack* otrack)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelSeltzerBerger<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack);
  cudaDeviceSynchronize();

  elapsedTime = timer.Stop();

  return elapsedTime;
}

} // end namespace vecphys
