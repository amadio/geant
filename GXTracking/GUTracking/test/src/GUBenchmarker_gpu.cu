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
    //[low:high] default energy limit = [10keV:1TeV]
    GUAliasSampler sampler(devStates,tid,1.e-2,1.e+6,100,200,table[kKleinNishina]);
    model.SetSampler(&sampler);
  }

  //  sampler.PrintTable();

  while (tid < nTrackSize) {
    model.Interact<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__
void KernelHybridCompton(Random_t* devStates,
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
    //[low:high] default energy limit = [10keV:1TeV]
    GUAliasSampler sampler(devStates,tid,1.e-2,1.e+6,100,200,table[kHybridKleinNishina]);
    model.SetSampler(&sampler);
  }

  while (tid < nTrackSize) {
    model.ModelInteract<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
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
    //[low:high] default energy limit = [2*electron_mass_c2:1TeV]
    GUAliasSampler sampler(devStates,tid,2.*electron_mass_c2,1.e+6,100,100,table[kBetheHeitler]);
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
    //[low:high] default energy limit = [0.1keV:1TeV] - should cut at Egamma/electron_mass_c2 > 50
    GUAliasSampler sampler(devStates,tid,1.e-4,1.e+6,100,200,table[kSauterGavrila]);
    model.SetSampler(&sampler);
  }

  while (tid < nTrackSize) {
    model.ModelInteract<kCuda>(itrack[tid],targetElements[tid],otrack[tid]);
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
    //[low:high] default energy limit = [0.1keV:1TeV]
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
    //[low:high] default energy limit = [0.1keV:1TeV]
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

void CudaKleinNishina(int blocksPerGrid, 
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
  vecphys::cuda::KernelKleinNishina<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
}

void CudaHybridCompton(int blocksPerGrid, 
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
  vecphys::cuda::KernelHybridCompton<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
}

void CudaBetheHeitler(int blocksPerGrid, 
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
  vecphys::cuda::KernelBetheHeitler<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
}

void CudaSauterGavrila(int blocksPerGrid, 
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
  vecphys::cuda::KernelSauterGavrila<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
}

void CudaMollerBhabha(int blocksPerGrid, 
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
  vecphys::cuda::KernelMollerBhabha<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
}

void CudaSeltzerBerger(int blocksPerGrid, 
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
  vecphys::cuda::KernelSeltzerBerger<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,sbData,nTrackSize,
    				       itrack, targetElements,otrack,sampleType);
}

} // end namespace vecphys
