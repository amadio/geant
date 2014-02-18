class GPFieldMap;
class GPPhysicsTable;
class GPVGeometry;
class GXFieldMap;
class GPPhysicsTable;
struct GXTrack;

#include <cuda.h>
#include <curand.h>
#include "random_kernel.h"

void tracking_gpu(curandState* devStates,
                  GPGeomManager *geomManager,
                  GXFieldMap *magMap,
                  GXTrack *track,
                  int *logVolumeIndices,
                  int *physVolumeIndices,
                  GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
                  size_t nTrackSize,
                  int NBLOCKS,
                  int NTHREADS,
                  cudaStream_t stream);

