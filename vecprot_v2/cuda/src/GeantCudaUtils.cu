// See also GeantCudaUtils.cxx

#include "GeantCudaUtils.h"

namespace Geant {
inline namespace cuda {

void CoprocessorBrokerInitConstant() {
   GEANT_CUDA_ERROR( cudaMemcpyToSymbol(gPropagator_fBmag, &(gPropagator->fBmag), sizeof(double), size_t(0), cudaMemcpyHostToDevice) );

   double tolerance = TGeoShape::Tolerance();
   GEANT_CUDA_ERROR( cudaMemcpyToSymbol(device_constant::gTolerance, &(tolerance), sizeof(double), size_t(0), cudaMemcpyHostToDevice) );
}

} // cuda
} // Geant
