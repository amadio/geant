#include "GeantFwd.h"

#include <cuda.h>
#include <curand.h>

namespace geant {
#ifdef VECCORE_CUDA
inline
#endif
namespace cuda {

class GeantTrack_v;
class GeantTaskData;

#ifdef VECCORE_CUDA
template <typename DataType, typename... ArgsTypes>
__global__ void MakeInstanceAtOnGpu(DataType *addr, ArgsTypes... params)
{
  DataType::MakeInstanceAt(addr, params...);
}

template <typename DataType, typename... ArgsTypes>
void MakeInstanceAt(DataType *addr, ArgsTypes... params)
{
  MakeInstanceAtOnGpu<<<1, 1>>>(addr, params...);
}

template <typename DataType, typename... ArgsTypes>
__global__ void MakeInstanceArrayAtOnGpu(DataType *addr, size_t nElements, size_t sizeOf, ArgsTypes... params)
{
  // nElements; number of element in the array.
  // sizeOf: actual size of each element.

  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  unsigned int idx = tid;
  while (idx < nElements) {
    // printf("Calling MakeInstanceAt from tid=%d for idx=%d\n",tid,idx);
    DataType::MakeInstanceAt(((char *)addr) + idx * sizeOf, params...);
    idx += blockDim.x * gridDim.x;
  }
}

template <typename DataType, typename... ArgsTypes>
void MakeInstanceArrayAt(DataType *addr, size_t nElements, size_t sizeOf, ArgsTypes... params)
{
  MakeInstanceArrayAtOnGpu<<<10, 1>>>(addr, nElements, sizeOf, params...);
}

#else
template <typename DataType, typename... ArgsTypes>
void MakeInstanceAt(DataType *addr, ArgsTypes... params);
template <typename DataType, typename... ArgsTypes>
void MakeInstanceArrayAt(DataType *addr, size_t nElements, size_t sizeOf, ArgsTypes... params);
#endif

} // cuda
} // Geant

int PropagateGeantTrack_gpu(vecgeom::cxx::DevicePtr<geant::cuda::GeantTaskData> &workSpace, size_t workspaceSizeOf,
                            size_t ntracks, vecgeom::cxx::DevicePtr<geant::cuda::GeantTrack_v> &input,
                            vecgeom::cxx::DevicePtr<geant::cuda::GeantTrack_v> &output,

                            int nBlocks, int nThreads, cudaStream_t stream);
namespace geant {
#ifdef VECCORE_CUDA
inline
#endif
namespace cuda {

int Clear_gpu(vecgeom::cxx::DevicePtr<geant::cuda::GeantTrack_v> &tracks, int nBlocks, int nThreads,
              cudaStream_t stream);
} // cuda
} // Geant
