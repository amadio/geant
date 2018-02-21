// See also GeantCudaUtils.cxx

#include "GeantCudaUtils.h"
#include "GeantPropagator.h"
#include "GeantTrack.h"
#include "GeantConfig.h"
#include "globals.h"

namespace Geant {
inline namespace cuda {

void CoprocessorBrokerInitConstant()
{
  double tolerance = 1e-7;
  GEANT_CUDA_ERROR(
      cudaMemcpyToSymbol(device_constant::gTolerance, &(tolerance), sizeof(double), size_t(0), cudaMemcpyHostToDevice));
}

} // cuda
} // Geant
