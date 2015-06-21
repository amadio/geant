//===------------------ Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantCudaUtils.h
 * @brief Utility routines for NVCC/CUDA compilation and setup
 * @author Philippe Canal
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_CUDAUTILS_H
#define GEANT_CUDAUTILS_H

#ifndef GEANT_CONFIG_H
#include "Geant/Config.h"
#endif

#include <cuda.h>
// #include "driver_types.h" // Required for cudaError_t type
#include "cuda_runtime.h" // Required for cudaGetErrorString

// This should be part of a global (configure time generated) header.
#ifndef VECGEOM_CUDA
#define VECGEOM_CUDA
#endif

#define GEANT_CUDA_ERROR( err ) (Geant::HandleCudaError( err, __FILE__, __LINE__ ))

namespace Geant {
   enum class EMsgLevel {
      kUnset    =  -1,
      kPrint    =   0,
      kInfo     =   1000,
      kWarning  =   2000,
      kError    =   3000,
      kBreak    =   4000,
      kSysError =   5000,
      kFatal    =   6000
   };
inline namespace GEANT_IMPL_NAMESPACE {

#ifndef GEANT_NVCC
   void ErrorHandlerImpl(EMsgLevel level, const char *location, const char *msgfmt, ...);
#endif

   template <typename... ArgsTypes>
   GEANT_CUDA_BOTH_CODE
   void MessageHandler(EMsgLevel level, const char *location, const char *msgfmt, ArgsTypes... params) {
#ifdef GEANT_NVCC
      const char *type = nullptr;
      switch(level) {
         case EMsgLevel::kPrint: type = "Print"; break;
         case EMsgLevel::kInfo: type = "Info"; break;
         case EMsgLevel::kWarning: type = "Warning"; break;
         case EMsgLevel::kError: type = "Error"; break;
         case EMsgLevel::kBreak: type = "Break"; break;
         case EMsgLevel::kSysError: type = "SysError"; break;
         case EMsgLevel::kFatal: type = "Fatal"; break;
         default: type = "Unknown Level"; break;
      }
      if (level == EMsgLevel::kPrint)
         printf("%s:",location);
      else
         printf("%s in <%s>:",type, location);
      printf(msgfmt,params...);
      printf("\n");
      if (level >= EMsgLevel::kFatal) {
#ifdef GEANT_CUDA_DEVICE_BUILD
         // Did not find a way to halt a kernel from within yet.
         //cudaDeviceReset();
         //cudaThreadExit();
         //throw("Fatal error in CUDA kernel");
#else
         exit( EXIT_FAILURE );
#endif
      }
#else
      // Currently we use the ROOT message handler on the host/gcc code.
      Geant::cxx::ErrorHandlerImpl(level,location,msgfmt,params...);
#endif
   }

   template <typename... ArgsTypes>
   GEANT_CUDA_BOTH_CODE
   void Error(const char *location, const char *msgfmt, ArgsTypes... params)
   {
      MessageHandler(EMsgLevel::kError,location,msgfmt, params...);
   }

   template <typename... ArgsTypes>
   GEANT_CUDA_BOTH_CODE
   void Fatal(const char *location, const char *msgfmt, ArgsTypes... params)
   {
      MessageHandler(EMsgLevel::kFatal,location,msgfmt, params...);
   }

   inline void HandleCudaError( cudaError_t err,
                                const char *file,
                                int line ) {
      if (err != cudaSuccess) {
         ::Geant::Fatal("Cuda","%s (%d) in %s at line %d\n", cudaGetErrorString( err ), err,
                        file, line );
      }
   }

} // GEANT_IMPL_NAMESPACE
} // Geant


#endif // GEANT_CUDAUTILS_H
