#ifndef GUBENCHMARKER_GPU_H
#define GUBENCHMARKER_GPU_H 1

#include "GUTrack.h"
#include "base/Global.h"
#include "GUAliasTableManager.h"
#include "Physics2DVector.h"
#include "SamplingMethod.h"

namespace vecphys {

// CUDA

void CudaKleinNishina(int blocksPerGrid,
		      int threadsPerBlock,
		      Random_t* devStates,
		      GUAliasTableManager** table,
		      Physics2DVector* sbData,
		      int nTrackSize,
		      GUTrack* itrack,
		      int* targetElements,
		      GUTrack* otrack,
		      SamplingMethod sampleType);

void CudaHybridCompton(int blocksPerGrid,
                       int threadsPerBlock,
		       Random_t* devStates,
		       GUAliasTableManager** table,
		       Physics2DVector* sbData,
		       int nTrackSize,
		       GUTrack* itrack,
		       int* targetElements,
		       GUTrack* otrack,
		       SamplingMethod sampleType);

void CudaBetheHeitler(int blocksPerGrid,
		      int threadsPerBlock,
		      Random_t* devStates,
		      GUAliasTableManager** table,
		      Physics2DVector* sbData,
		      int nTrackSize,
		      GUTrack* itrack,
		      int* targetElements,
		      GUTrack* otrack,
		      SamplingMethod sampleType);

void CudaSauterGavrila(int blocksPerGrid,
		       int threadsPerBlock,
		       Random_t* devStates,
		       GUAliasTableManager** table,
		       Physics2DVector* sbData,
		       int nTrackSize,
		       GUTrack* itrack,
		       int* targetElements,
		       GUTrack* otrack,
		       SamplingMethod sampleType);

void CudaMollerBhabha(int blocksPerGrid,
		      int threadsPerBlock,
		      Random_t* devStates,
		      GUAliasTableManager** table,
		      Physics2DVector* sbData,
		      int nTrackSize,
		      GUTrack* itrack,
		      int* targetElements,
		      GUTrack* otrack,
		      SamplingMethod sampleType);

void CudaSeltzerBerger(int blocksPerGrid,
		       int threadsPerBlock,
		       Random_t* devStates,
		       GUAliasTableManager** table,
		       Physics2DVector* sbData,
		       int nTrackSize,
		       GUTrack* itrack,
		       int* targetElements,
		       GUTrack* otrack,
		       SamplingMethod sampleType);

typedef void (*CudaKernelFunc_t)(int blocksPerGrid,
				 int threadsPerBlock,
				 Random_t* devStates,
				 GUAliasTableManager** table,
				 Physics2DVector* sbData,
				 int nTrackSize,
				 GUTrack* itrack,
				 int* targetElements,
				 GUTrack* otrack,
				 SamplingMethod sampleType);

CudaKernelFunc_t CudaKernelFunc[] = {CudaKleinNishina,
                                     CudaHybridCompton,
                                     CudaBetheHeitler,
                                     CudaSauterGavrila,
                                     CudaMollerBhabha,
                                     CudaSeltzerBerger};

} // end namespace vecphys

#endif
