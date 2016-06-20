#ifndef PROCESSBENCHMARKER_GPU_H
#define PROCESSBENCHMARKER_GPU_H 1

#include "base/VecPhys.h"
#include "GUTrack.h"
#include "EmProcess.h"

namespace vecphys {

// CUDA
void CudaPhotonProcess(int blocksPerGrid, int threadsPerBlock, Random_t *devStates, CrossSectionData *table,
                       int nTrackSize, GUTrack *itrack, int *targetElements);

void CudaElectronProcess(int blocksPerGrid, int threadsPerBlock, Random_t *devStates, CrossSectionData *table,
                         int nTrackSize, GUTrack *itrack, int *targetElements);

typedef void (*CudaKernelFunc_t)(int blocksPerGrid, int threadsPerBlock, Random_t *devStates, CrossSectionData *table,
				 int nTrackSize, GUTrack *itrack, int *targetElements);

CudaKernelFunc_t CudaKernelFunc[] = {CudaPhotonProcess, CudaElectronProcess};

} // end namespace vecphys

#endif
