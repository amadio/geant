class GPFieldMap;
class GPPhysicsTable;
class GPVGeometry;
class GXFieldMap;
class GPPhysicsTable;
struct GXTrack;
class GXTrackLiason;
class GPGeomManager;
class PPhysics2DVector;
class GPPhysics2DVector;

class GeantTaskData;

#include <cuda.h>
#include <curand.h>
#include "random_kernel.h"

namespace Geant {
namespace cuda {

#ifdef GEANT_NVCC
   template <typename DataType, typename... ArgsTypes>
   __global__
   void MakeInstanceAtOnGpu(DataType *addr, ArgsTypes... params) {
      DataType::MakeInstanceAt(addr, params...);
   }

   template <typename DataType, typename... ArgsTypes>
   void MakeInstanceAt(DataType *addr, ArgsTypes... params)
   {
      MakeInstanceAtOnGpu<<<1, 1>>>(addr, params...);
   }
#else
   template <typename DataType, typename... ArgsTypes>
      void MakeInstanceAt(DataType *addr, ArgsTypes... params);
#endif

} // cuda
} // Geant


int tracking_gpu(curandState* devStates,
                 size_t nSteps,
                 size_t nElectrons,
                 GXTrack *track, GXTrack * /* altTrack */,
                 int *logVolumeIndices,
                 int *physVolumeIndices,
                 GXTrack *secondaries, int *secStackSize,

                 int *scratch,
                 GXTrackLiason *trackScratch,

                 GPGeomManager *geomManager,
                 GXFieldMap *magMap,
                 GPPhysicsTable *physicsTable,
                 GPPhysics2DVector *seltzerBergerTable,

                 int nBlocks, int nThreads,
                 cudaStream_t stream);


int electron_gpu(curandState* devStates,
                 size_t nSteps,
                 size_t nElectrons,
                 GXTrack *track, GXTrack * /* altTrack */,
                 int *logVolumeIndices,
                 int *physVolumeIndices,
                 GXTrack *secondaries, int *secStackSize,

                 int *scratch,
                 GXTrackLiason *trackScratch,

                 GPGeomManager *geomManager,
                 GXFieldMap *magMap,
                 GPPhysicsTable *physicsTable,
                 GPPhysics2DVector *seltzerBergerTable,

                 int nBlocks, int nThreads,
                 cudaStream_t stream);

int electron_multistage_gpu(curandState* devStates,
                 size_t nSteps,
                 size_t nElectrons,
                 GXTrack *track, GXTrack * /* altTrack */,
                 int *logVolumeIndices,
                 int *physVolumeIndices,
                 GXTrack *secondaries, int *secStackSize,

                 int *scratch,
                 GXTrackLiason *trackScratch,

                 GPGeomManager *geomManager,
                 GXFieldMap *magMap,
                 GPPhysicsTable *physicsTable,
                 GPPhysics2DVector *seltzerBergerTable,

                 int nBlocks, int nThreads,
                 cudaStream_t stream);

int PropagateGeantTrack_gpu(vecgeom::cxx::DevicePtr<GeantTaskData> &workSpace,
                            size_t ntracks,
                            vecgeom::cxx::DevicePtr<GeantTrack_v> &input,
                            vecgeom::cxx::DevicePtr<GeantTrack_v> &output,

                            int nBlocks, int nThreads,
                            cudaStream_t stream);
