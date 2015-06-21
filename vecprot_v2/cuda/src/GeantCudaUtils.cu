// See also GeantCudaUtils.cxx

#include "GeantCudaUtils.h"

namespace Geant {
inline namespace cuda {

void Fatal(const char *location, const char *va_(fmt), ...)
{
   // Use this function in case of a fatal error. It will abort the program.

   va_list ap;
   va_start(ap,va_(fmt));
   printf("Fatal in <%s>:",location);
   printf(va_(fmt), ap);
   printf("\n");
   va_end(ap);
#ifdef GEANT_CUDA_DEVICE_BUILD
   //cudaDeviceReset();
   cudaThreadExit();
#else
   exit( EXIT_FAILURE );
#endif
}

void CoprocessorBrokerInitConstant() {
   GEANT_CUDA_ERROR( cudaMemcpyToSymbol(gPropagator_fBmag, &(gPropagator->fBmag), sizeof(double), size_t(0), cudaMemcpyHostToDevice) );

   double tolerance = TGeoShape::Tolerance();
   GEANT_CUDA_ERROR( cudaMemcpyToSymbol(device_constant::gTolerance, &(tolerance), sizeof(double), size_t(0), cudaMemcpyHostToDevice) );
}

} // cuda
} // Geant
