#include "Geant/Error.h"
#ifdef USE_ROOT
#include "TError.h"
#endif
#include <stdarg.h>

namespace Geant {
inline namespace cxx {

// Code to be compiled only by gcc (i.e. not nvcc).

#ifdef USE_ROOT
void ErrorHandlerImpl(EMsgLevel level, const char *location, const char *fmt, ...)
{
    // Currently we use the ROOT message handler on the host/gcc code.

   va_list ap;
   va_start(ap,fmt);
   ::ErrorHandler((int)level,location, fmt, ap);
   va_end(ap);
}
#else
void ErrorHandlerImpl(EMsgLevel level, const char *location, const char *msgfmt, ...)
{
   // Trivial implementation
   va_list args;
   va_start(args, msgfmt);

   if (level > EMsgLevel::kPrint || (location==nullptr || location[0]=='\0')) {
     fprintf(stdout,"Geant Message level %d at %s:",(int)level,location);
   }
   vfprintf(stdout, msgfmt, args);
   fprintf(stdout, "\n");
   va_end(args);
}
#endif

} // cxx
} // Geant
