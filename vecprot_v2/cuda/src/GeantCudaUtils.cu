// See also GeantCudaUtils.cxx

#include "GeantCudaUtils.h"
#include "GeantPropagator.h"
#include "GeantTrack.h"
#include "GeantConfig.h"
#ifndef USE_VECGEOM_NAVIGATOR
#include "TGeoShape.h"
#endif
#include "globals.h"

namespace Geant {
inline namespace cuda {

void CoprocessorBrokerInitConstant()
{
#ifndef USE_VECGEOM_NAVIGATOR
  double tolerance = TGeoShape::Tolerance();
#else
  double tolerance = 1e-7;
#endif
  GEANT_CUDA_ERROR(
      cudaMemcpyToSymbol(device_constant::gTolerance, &(tolerance), sizeof(double), size_t(0), cudaMemcpyHostToDevice));
}

} // cuda
} // Geant
