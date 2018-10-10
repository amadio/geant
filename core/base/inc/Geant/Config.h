/// \file Geant/Config.h

#ifndef GEANT_CONFIG_H
#define GEANT_CONFIG_H

#if __cplusplus < 201103L && !defined(__NVCC__)
#error "GeantV requires C++11"
#endif

// Include global definitions from VecCore
#include "base/Global.h"

// Inlining
#ifdef __INTEL_COMPILER
#define GEANT_FORCE_INLINE inline
#else
#if (defined(__GNUC__) || defined(__GNUG__)) && !defined(__clang__) && !defined(__NO_INLINE__) && \
    !defined(GEANT_NOINLINE)
#define GEANT_FORCE_INLINE inline __attribute__((always_inline))
#else
// Clang or forced inlining is disabled ( by falling back to compiler decision )
#define GEANT_FORCE_INLINE inline
#endif
#endif

// Utility macros
#define GEANT_BIT(n)       (1ULL << (n))
#define GEANT_SETBIT(n,i)  ((n) |= BIT(i))
#define GEANT_CLRBIT(n,i)  ((n) &= ~BIT(i))
#define GEANT_TESTBIT(n,i) ((bool)(((n) & BIT(i)) != 0))

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
// to geant::cxx::gTolerance and geant::cuda::gTolerance respectively.
// In nvcc device code, geant::cuda::gTolerance is aliased to
// device_constant::gTolerance.

#ifndef VECCORE_CUDA

#define GEANT_IMPL_NAMESPACE cxx

#define GEANT_DECLARE_CONSTANT(type, name) \
  namespace host_constant {                \
  extern const type name;                  \
  }                                        \
  using host_constant::name
#else

#define GEANT_IMPL_NAMESPACE cuda

#ifdef VECCORE_CUDA_DEVICE_COMPILATION
#ifdef CUDA_SEP_COMP
#define GEANT_DECLARE_CONSTANT(type, name) \
  namespace host_constant {                \
  extern const type name;                  \
  }                                        \
  namespace device_constant {              \
  extern __constant__ type name;           \
  }                                        \
  using device_constant::name
#else // CUDA_SEP_COMP
#define GEANT_DECLARE_CONSTANT(type, name) \
  namespace host_constant {                \
  extern const type name;                  \
  }                                        \
  namespace device_constant {              \
  __constant__ type name;                  \
  }                                        \
  using device_constant::name
#endif // CUDA_SEP_COMP
#else  // VECCORE_CUDA_DEVICE_COMPILATION
#ifdef CUDA_SEP_COMP
#define GEANT_DECLARE_CONSTANT(type, name) \
  namespace host_constant {                \
  extern const type name;                  \
  }                                        \
  namespace device_constant {              \
  extern __constant__ type name;           \
  }                                        \
  using host_constant::name
#else // CUDA_SEP_COMP
#define GEANT_DECLARE_CONSTANT(type, name) \
  namespace host_constant {                \
  extern const type name;                  \
  }                                        \
  namespace device_constant {              \
  __constant__ type name;                  \
  }                                        \
  using host_constant::name
#endif // CUDA_SEP_COMP
#endif // Device build or not
#endif // gcc or nvcc

#ifndef VECCORE_CUDA

#define GEANT_DEVICE_DECLARE_CONV(NS, classOrStruct, X) \
  namespace NS {                                        \
  namespace cuda {                                      \
  classOrStruct X;                                      \
  }                                                     \
  inline namespace cxx {                                \
  classOrStruct X;                                      \
  }                                                     \
  }                                                     \
  namespace vecgeom {                                   \
  template <>                                           \
  struct kCudaType<NS::cxx::X> {                        \
    using type_t = NS::cuda::X;                         \
  };                                                    \
  }                                                     \
  class __QuietSemi

#else

#define GEANT_DEVICE_DECLARE_CONV(NS, classOrStruct, X) class __QuietSemi

#endif

#endif
