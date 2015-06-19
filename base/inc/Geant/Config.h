/// \file Geant/Config.h

#ifndef GEANT_CONFIG_H
#define GEANT_CONFIG_H

#if __cplusplus < 201103L && !defined(__NVCC__)
#error "GeantV requires C++11"
#endif

#if (defined(__CUDACC__) || defined(__NVCC__))

  // Compiling with nvcc
  #define GEANT_NVCC
  #ifdef __CUDA_ARCH__
     #define GEANT_CUDA_DEVICE_BUILD
  #endif
  #define GEANT_CUDA_DEVICE_CODE __device__
  #define GEANT_CUDA_HOST_CODE __host__
  #define GEANT_CUDA_BOTH_CODE __host__ __device__
  #define GEANT_DEVICE_CONSTANT __constant__

  #define GEANT_IMPL_NAMESPACE cuda

#else
  // Not compiling with NVCC
   #define GEANT_CUDA_DEVICE_CODE
   #define GEANT_CUDA_HOST_CODE
   #define GEANT_CUDA_BOTH_CODE 
   #define GEANT_DEVICE_CONSTANT extern const

   #define GEANT_IMPL_NAMESPACE cxx

#endif

#endif
