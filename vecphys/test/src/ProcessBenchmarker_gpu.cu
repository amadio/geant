#include "PhotonProcess.h"
#include "ElectronProcess.h"

namespace vecphys {
inline namespace cuda {

__global__
void KernelPhotonProcess(Random_t* devStates,
                         CrossSectionData* table,
			 int nTrackSize, 
                         GUTrack* itrack, 
			 int* materialIndex) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  PhotonProcess process(devStates,tid,table);

  while (tid < nTrackSize) {
    process.GetStepLengthAndProcess<backend::Scalar>(itrack[tid], materialIndex[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

__global__
void KernelElectronProcess(Random_t* devStates,
                           CrossSectionData* table,
		 	   int nTrackSize, 
                           GUTrack* itrack, 
			   int* materialIndex) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  ElectronProcess process(devStates,tid,table);

  while (tid < nTrackSize) {
    process.GetStepLengthAndProcess<backend::Scalar>(itrack[tid], materialIndex[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

} // end namespace cuda

// Cuda wrapper

void CudaPhotonProcess(int blocksPerGrid, 
		       int threadsPerBlock,
		       Random_t* devStates,
		       CrossSectionData* table,
 		       int nTrackSize, 
		       GUTrack* itrack, 
		       int* targetElements) 
{
  vecphys::cuda::KernelPhotonProcess<<<blocksPerGrid, threadsPerBlock>>>(
				       devStates,table,nTrackSize,itrack,targetElements);
}

void CudaElectronProcess(int blocksPerGrid, 
		         int threadsPerBlock,
		         Random_t* devStates,
		         CrossSectionData* table,
 		         int nTrackSize, 
		         GUTrack* itrack, 
		         int* targetElements) 
{
  vecphys::cuda::KernelElectronProcess<<<blocksPerGrid, threadsPerBlock>>>(
				         devStates,table,nTrackSize,itrack,targetElements);
}

} // end namespace vecphys
