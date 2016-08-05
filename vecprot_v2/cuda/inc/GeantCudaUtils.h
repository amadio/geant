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

#ifndef GEANT_ERROR_H
#include "Geant/Error.h"
#endif

#include <cuda.h>
// #include "driver_types.h" // Required for cudaError_t type
#include "cuda_runtime.h" // Required for cudaGetErrorString

#define GEANT_CUDA_ERROR(err) (Geant::HandleCudaError(err, __FILE__, __LINE__))

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

inline void HandleCudaError(cudaError_t err, const char *file, int line)
{
  if (err != cudaSuccess) {
    ::Geant::Fatal("Cuda", "%s (%d) in %s at line %d\n", cudaGetErrorString(err), err, file, line);
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANT_CUDAUTILS_H
