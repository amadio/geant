#include "Geant/Error.h"
#include "TError.h"
#include <stdarg.h>
#include "Varargs.h"

namespace Geant {
inline namespace cxx {

// Code to be compiled only by gcc (i.e. not nvcc).

void ErrorHandlerImpl(EMsgLevel level, const char *location, const char *va_(fmt), ...)
{
    // Currently we use the ROOT message handler on the host/gcc code.

   va_list ap;
   va_start(ap,va_(fmt));
   ::ErrorHandler((int)level,location, va_(fmt), ap);
   va_end(ap);
}

} // cxx
} // Geant
