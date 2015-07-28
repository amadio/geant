#ifndef GUBENCHMARKER_GPU_H
#define GUBENCHMARKER_GPU_H 1

#include "GUTrack.h"
#include "base/Global.h"
#include "GUAliasTableManager.h"
#include "Physics2DVector.h"

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
			   GUTrack* otrack);

Precision CudaBetheHeitler(int blocksPerGrid, 
                           int threadsPerBlock,
                           Random_t* devStates,
                           GUAliasTableManager** table,
                           Physics2DVector* sbData,
			   int nTrackSize, 
                           GUTrack* itrack, 
			   int* targetElements, 
			   GUTrack* otrack);

Precision CudaSauterGavrila(int blocksPerGrid, 
                            int threadsPerBlock,
                            Random_t* devStates,
                            GUAliasTableManager** table,
                            Physics2DVector* sbData,
			    int nTrackSize, 
                            GUTrack* itrack, 
			    int* targetElements, 
			    GUTrack* otrack);

Precision CudaMollerBhabha(int blocksPerGrid, 
                           int threadsPerBlock,
                           Random_t* devStates,
                           GUAliasTableManager** table,
                           Physics2DVector* sbData,
			   int nTrackSize, 
                           GUTrack* itrack, 
			   int* targetElements, 
			   GUTrack* otrack);

Precision CudaSeltzerBerger(int blocksPerGrid, 
                            int threadsPerBlock,
                            Random_t* devStates,
                            GUAliasTableManager** table,
                            Physics2DVector* sbData,
			    int nTrackSize, 
                            GUTrack* itrack, 
			    int* targetElements, 
			    GUTrack* otrack);

typedef Precision (*CudaKernelFunc_t)(int blocksPerGrid, 
                                      int threadsPerBlock,
                                      Random_t* devStates,
                                      GUAliasTableManager** table,
                                      Physics2DVector* sbData,
			              int nTrackSize, 
                                      GUTrack* itrack, 
			              int* targetElements, 
			              GUTrack* otrack);

CudaKernelFunc_t CudaKernelFunc[] = {CudaKleinNishina,
                                     CudaBetheHeitler, 
                                     CudaSauterGavrila, 
                                     CudaMollerBhabha, 
                                     CudaSeltzerBerger}; 

} // end namespace vecphys

#endif
