// See also GeantCudaUtils.cxx

#include "GeantCudaUtils.h"
#include "Propagator.h"
#include "Geant/Track.h"
#include "GeantConfig.h"

namespace geant {
inline namespace cuda {

void CoprocessorBrokerInitConstant()
{
  double tolerance = 1e-7;
  GEANT_CUDA_ERROR(
      cudaMemcpyToSymbol(device_constant::gTolerance, &(tolerance), sizeof(double), size_t(0), cudaMemcpyHostToDevice));
}

} // cuda
} // Geant
