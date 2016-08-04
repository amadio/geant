// See also GeantCudaUtils.cu

#include "GeantCudaUtils.h"

#ifdef USE_ROOT
#include "TError.h"
#include "Varargs.h"
#endif
#include <stdarg.h>

namespace Geant {
inline namespace cxx {

// Code to be compiled only by gcc (i.e. not nvcc).

} // cxx
} // Geant
