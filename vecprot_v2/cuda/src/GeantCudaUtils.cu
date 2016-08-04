// See also GeantCudaUtils.cxx

#include "GeantCudaUtils.h"
#include "GeantPropagator.h"
#include "GeantTrack.h"
#ifdef USE_ROOT
#include "TGeoShape.h"
#endif
#include "globals.h"

namespace Geant {
inline namespace cuda {

void CoprocessorBrokerInitConstant()
{
  GEANT_CUDA_ERROR(
      cudaMemcpyToSymbol(gPropagator_fBmag, &(gPropagator->fBmag), sizeof(double), size_t(0), cudaMemcpyHostToDevice));

#ifdef USE_ROOT
  double tolerance = TGeoShape::Tolerance();
#else
  double tolerance = 1e-7;
#endif
  GEANT_CUDA_ERROR(
      cudaMemcpyToSymbol(device_constant::gTolerance, &(tolerance), sizeof(double), size_t(0), cudaMemcpyHostToDevice));
}

} // cuda
} // Geant
