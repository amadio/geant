#ifndef GUBENCHMARKER_GPU_H
#define GUBENCHMARKER_GPU_H 1

#include "GUTrack.h"
#include "base/Global.h"
#include "GUAliasTableManager.h"
#include "Physics2DVector.h"
#include "SamplingMethod.h"

namespace vecphys {

// CUDA

Precision CudaKleinNishina(int blocksPerGrid, 
                           int threadsPerBlock,
                           Random_t* devStates,
                           GUAliasTableManager** table,
                           Physics2DVector* sbData,
			   int nTrackSize, 
                           GUTrack* itrack, 
			   int* targetElements, 
			   GUTrack* otrack,
                           SamplingMethod sampleType);

Precision CudaBetheHeitler(int blocksPerGrid, 
                           int threadsPerBlock,
                           Random_t* devStates,
                           GUAliasTableManager** table,
                           Physics2DVector* sbData,
			   int nTrackSize, 
                           GUTrack* itrack, 
			   int* targetElements, 
			   GUTrack* otrack,
                           SamplingMethod sampleType);

Precision CudaSauterGavrila(int blocksPerGrid, 
                            int threadsPerBlock,
                            Random_t* devStates,
                            GUAliasTableManager** table,
                            Physics2DVector* sbData,
			    int nTrackSize, 
                            GUTrack* itrack, 
			    int* targetElements, 
			    GUTrack* otrack,
                            SamplingMethod sampleType);

Precision CudaMollerBhabha(int blocksPerGrid, 
                           int threadsPerBlock,
                           Random_t* devStates,
                           GUAliasTableManager** table,
                           Physics2DVector* sbData,
			   int nTrackSize, 
                           GUTrack* itrack, 
			   int* targetElements, 
			   GUTrack* otrack,
                           SamplingMethod sampleType);

Precision CudaSeltzerBerger(int blocksPerGrid, 
                            int threadsPerBlock,
                            Random_t* devStates,
                            GUAliasTableManager** table,
                            Physics2DVector* sbData,
			    int nTrackSize, 
                            GUTrack* itrack, 
			    int* targetElements, 
			    GUTrack* otrack,
			    SamplingMethod sampleType);

typedef Precision (*CudaKernelFunc_t)(int blocksPerGrid, 
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
                                     CudaBetheHeitler, 
                                     CudaSauterGavrila, 
                                     CudaMollerBhabha, 
                                     CudaSeltzerBerger}; 

} // end namespace vecphys

#endif
