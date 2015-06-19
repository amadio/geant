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
inline namespace GEANT_IMPL_NAMESPACE {

   void Fatal(const char *location, const char *msgfmt, ...)
   #if defined(__GNUC__)
      __attribute__((format(printf, 2, 3)))
   #endif
      ;

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
