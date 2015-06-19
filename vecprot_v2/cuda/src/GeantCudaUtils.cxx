// See also GeantCudaUtils.cu

#include "GeantCudaUtils.h"

#include "TError.h"
#include <stdarg.h>
#include "Varargs.h"

namespace Geant {
inline namespace cxx {

void Fatal(const char *location, const char *va_(fmt), ...)
{
   // Use this function in case of a fatal error. It will abort the program.

   va_list ap;
   va_start(ap,va_(fmt));
   ::Fatal(location, va_(fmt), ap);
   va_end(ap);
#ifdef GEANT_CUDA_DEVICE_BUILD
   cudaThreadExit();
#else
   exit( EXIT_FAILURE );
#endif
}

} // cxx
} // Geant
