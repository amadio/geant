#include "base/Stopwatch.h"

//#include "GVComptonKleinNishina.h"
#include "ComptonKleinNishina.h"
#include "ConversionBetheHeitler.h"
#include "PhotoElectronSauterGavrila.h"
#include "IonisationMoller.h"
#include "BremSeltzerBerger.h"
#include "SamplingMethod.h"

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
			GUTrack* otrack,
                        SamplingMethod sampleType)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  ComptonKleinNishina model(devStates,tid,0);
  model.SetSamplingMethod(sampleType);
  if(model.GetSamplingMethod() == SamplingMethod::kAlias) { 
    GUAliasSampler sampler(devStates,tid,1.e-4,1.e+6,100,100,table[kKleinNishina]);
    model.SetSampler(&sampler);
  }

  //  sampler.PrintTable();

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
			GUTrack* otrack,
                        SamplingMethod sampleType)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  ConversionBetheHeitler model(devStates,tid,0);
  model.SetSamplingMethod(sampleType);
  if(model.GetSamplingMethod() == SamplingMethod::kAlias) { 
    GUAliasSampler sampler(devStates,tid,1.e-4,1.e+6,100,100,table[kBetheHeitler]);
    model.SetSampler(&sampler);
  }

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
			 GUTrack* otrack,
                         SamplingMethod sampleType)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  PhotoElectronSauterGavrila model(devStates,tid,0);
  model.SetSamplingMethod(sampleType);
  if(model.GetSamplingMethod() == SamplingMethod::kAlias) { 
    GUAliasSampler sampler(devStates,tid,1.e-4,1.e+6,100,100,table[kSauterGavrila]);
    model.SetSampler(&sampler);
  }

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
			GUTrack* otrack,
                        SamplingMethod sampleType)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  IonisationMoller model(devStates,tid,0);
  model.SetSamplingMethod(sampleType);
  if(model.GetSamplingMethod() == SamplingMethod::kAlias) { 
    GUAliasSampler sampler(devStates,tid,1.e-4,1.e+6,100,100,table[kMollerBhabha]);
    model.SetSampler(&sampler);
  }

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
			 GUTrack* otrack,
                         SamplingMethod sampleType)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  BremSeltzerBerger model(devStates,tid,0,sbData);
  model.SetSamplingMethod(sampleType);
  if(model.GetSamplingMethod() == SamplingMethod::kAlias) { 
    GUAliasSampler sampler(devStates,tid,1.e-4,1.e+6,100,100,table[kSeltzerBerger]);
    model.SetSampler(&sampler);
  }

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
			   GUTrack* otrack,
                           SamplingMethod sampleType)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelKleinNishina<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
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
			   GUTrack* otrack,
                           SamplingMethod sampleType)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelBetheHeitler<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
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
			    GUTrack* otrack,
                            SamplingMethod sampleType)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelSauterGavrila<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
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
			   GUTrack* otrack,
                           SamplingMethod sampleType)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelMollerBhabha<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
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
			    GUTrack* otrack,
                            SamplingMethod sampleType)
{
  static Stopwatch timer;
  Precision elapsedTime = 0.0;

  timer.Start();

  vecphys::cuda::KernelSeltzerBerger<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
  cudaDeviceSynchronize();

  elapsedTime = timer.Stop();

  return elapsedTime;
}

} // end namespace vecphys
