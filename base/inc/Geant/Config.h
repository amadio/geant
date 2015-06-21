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

//////////////////////////////////////////
// Declaration for constant define in the
// device constant section. Use:
//    GEANT_DECLARE_CONST(double,gTolerance);
//
// This will declare the following:
//    extern double host_constant::gTolerance;
// and only in nvcc
//    __constant device_constant::gTolerance;
//
// In gcc and nvcc host code host_constant::gTolerance is aliased
// to Geant::cxx::gTolerance and Geant::cuda::gTolerance respectively.
// In nvcc device code, Geant::cuda::gTolerance is aliased to
// device_constant::gTolerance.

#ifndef GEANT_NVCC
#define GEANT_DECLARE_CONSTANT(type,name) \
   namespace host_constant { \
      extern const type name; \
   } \
   using host_constant::name
#else
#ifdef GEANT_CUDA_DEVICE_BUILD
#define GEANT_DECLARE_CONSTANT(type,name) \
   namespace host_constant { \
      extern const type gTolerance; \
   } \
   namespace device_constant { \
      __constant__ type name; \
   } \
   using device_constant::gTolerance
#else
#define GEANT_DECLARE_CONSTANT(type,name) \
   namespace host_constant { \
      extern const type name; \
   } \
   namespace device_constant { \
      __constant__ type name; \
   } \
   using host_constant::gTolerance
#endif // Device build or not
#endif // gcc or nvcc

#endif
