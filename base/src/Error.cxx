#include "Geant/Error.h"
#ifdef USE_ROOT
#include "TError.h"
#include "Varargs.h"
#endif
#include <stdarg.h>

namespace Geant {
inline namespace cxx {

// Code to be compiled only by gcc (i.e. not nvcc).

#ifdef USE_ROOT
void ErrorHandlerImpl(EMsgLevel level, const char *location, const char *va_(fmt), ...)
{
    // Currently we use the ROOT message handler on the host/gcc code.

   va_list ap;
   va_start(ap,va_(fmt));
   ::ErrorHandler((int)level,location, va_(fmt), ap);
   va_end(ap);
}
#else
void ErrorHandlerImpl(EMsgLevel level, const char *location, const char *va_(msgfmt), ...)
{
   // Trivial implementation
   va_list args;
   va_start(args, va_(msgfmt));

   if (level > EMsgLevel::kPrint || (location==nullptr || location[0]=='\0')) {
    vfprintf(stdout,"Geant Message level %d at %s:",level,location);
   }
   vfprintf(stdout, msgfmt, args);
   vfprintf(stdout, "\n");
   va_end(ap);
}
#endif

} // cxx
} // Geant
