/// \file Global.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)
//
// simple version based on VecGeom - syjun
//
#ifndef VECPHYS_BASE_GLOBAL_H_
#define VECPHYS_BASE_GLOBAL_H_

#include <cassert>
#include <cmath>
#include <float.h>
#include <limits>
#include <stdio.h>
#include <memory>

#define VECPHYS

#if __cplusplus >= 201103L
  #define VECPHYS_STD_CXX11
#endif

#if (defined(__CUDACC__) || defined(__NVCC__))
  // Compiling with nvcc
  #define VECPHYS_NVCC
  #ifdef __CUDA_ARCH__
    // Compiling device code
    #define VECPHYS_NVCC_DEVICE
  #endif
  #define HAVECUDANAMESPACE
  #define VECPHYS_IMPL_NAMESPACE cuda
  #define VECPHYS_NAMESPACE ::vecphys
  #define VECPHYS_CUDA_HEADER_HOST __host__
  #define VECPHYS_CUDA_HEADER_DEVICE __device__
  #define VECPHYS_CUDA_HEADER_BOTH __host__ __device__
  #define VECPHYS_CUDA_HEADER_GLOBAL __global__
  #define VECPHYS_ALIGNED __align__((64))
  #define VECPHYS_HOST_FORWARD_DECLARE(X) namespace cxx { X }
  #define VECPHYS_DEVICE_FORWARD_DECLARE(X)
  #define VECPHYS_DEVICE_DECLARE_CONV(X)
  #define VECPHYS_DEVICE_DECLARE_CONV_TEMPLATE(X,ArgType)
  #undef VECPHYS_VC
  #undef VECPHYS_VC_ACCELERATION
  #undef VECPHYS_CILK
  #undef VECPHYS_ROOT
  #undef VECPHYS_USOLIDS
  #undef VECPHYS_GEANT4
  #undef VECPHYS_BENCHMARK
#else
  // Not compiling with NVCC
#define HAVENORMALNAMESPACE
  #define VECPHYS_IMPL_NAMESPACE cxx
  #define VECPHYS_NAMESPACE ::vecphys
  #define VECPHYS_CUDA_HEADER_HOST
  #define VECPHYS_CUDA_HEADER_DEVICE
  #define VECPHYS_CUDA_HEADER_BOTH
  #define VECPHYS_CUDA_HEADER_GLOBAL
  #ifdef VECPHYS_CUDA
    // CUDA is enabled, but currently compiling regular C++ code.
    // This enables methods that interface between C++ and CUDA environments
    #define VECPHYS_CUDA_INTERFACE
  #endif
  namespace vecphys {
     template <typename DataType> struct kCudaType;
     template <typename DataType> using CudaType_t = typename kCudaType<DataType>::type_t;
     template <> struct kCudaType<float> { using type_t = float; };
     template <> struct kCudaType<double> { using type_t = double; };
     template <> struct kCudaType<int> { using type_t = int; };
  }
  #define VECPHYS_HOST_FORWARD_DECLARE(X)
  #define VECPHYS_DEVICE_FORWARD_DECLARE(X)  namespace cuda { X }

  #define VECPHYS_DEVICE_DECLARE_CONV(X) \
     namespace cuda { class X; } \
     inline namespace cxx  { class X; } \
     template <> struct kCudaType<cxx::X> { using type_t = cuda::X; };

  #define VECPHYS_DEVICE_DECLARE_CONV_TEMPLATE(X,ArgType) \
     namespace cuda { template <ArgType Arg> class X; } \
     inline namespace cxx  { template <ArgType Arg> class X; } \
     template <ArgType Arg> struct kCudaType<cxx::X<Arg> > \
     { using type_t = cuda::X<CudaType_t<Arg> >; };

/* Instead of multiple macro, when we have auto expansion of Template pack we could use:
template <typename... Arguments>
struct kCudaType<cxx::BoxImplementation<Arguments...>  >
   { using type_t = typename cuda::BoxImplementation<CudaType_t<Arguments...> >; };
*/
#endif

#ifdef __INTEL_COMPILER
  // Compiling with icc
  #define VECPHYS_INTEL
  #define VECPHYS_INLINE inline
#else
  // Functionality of <mm_malloc.h> is automatically included in icc
  #include <mm_malloc.h>
  #if (defined(__GNUC__) || defined(__GNUG__)) && !defined(__clang__) && !defined(__NO_INLINE__) && !defined( VECPHYS_NOINLINE )
    #define VECPHYS_INLINE inline __attribute__((always_inline))
    #ifndef VECPHYS_NVCC
      #define VECPHYS_ALIGNED __attribute__((aligned(64)))
    #endif
  #else
#pragma message "forced inlining disabled"
  // Clang or forced inlining is disabled ( by falling back to compiler decision )
    #define VECPHYS_INLINE inline
    #ifndef VECPHYS_NVCC
      #define VECPHYS_ALIGNED
    #endif
  #endif
#endif

#ifndef NULL
  #define NULL nullptr
#endif

// Allow constexpr variables and functions if possible
#ifdef VECPHYS_STD_CXX11
  #define VECPHYS_CONSTEXPR constexpr
  #define VECPHYS_CONSTEXPR_RETURN constexpr
#else
  #define VECPHYS_CONSTEXPR const
  #define VECPHYS_CONSTEXPR_RETURN
#endif

// Qualifier(s) of global constants
#ifndef VECPHYS_NVCC
  #define VECPHYS_GLOBAL constexpr
#else
  #define VECPHYS_GLOBAL static __constant__ const
#endif

namespace vecphys {
#ifdef VECPHYS_FLOAT_PRECISION
typedef float Precision;
#else
typedef double Precision;
#endif

typedef int Inside_t;

}

#ifdef __CUDACC__
#include <cuda.h>
#include <curand_kernel.h>
typedef curandState Random_t;
#else
typedef int         Random_t;
#endif

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

#ifndef VECPHYS_NVCC
   using std::unique_ptr;

#else

   template <typename T>
   class unique_ptr {
      T *fValue;
   public:
     VECPHYS_CUDA_HEADER_BOTH
     unique_ptr(T *in) : fValue(in) {}

     VECPHYS_CUDA_HEADER_BOTH
     ~unique_ptr() { delete fValue; }

     VECPHYS_CUDA_HEADER_BOTH
     T* operator->() { return fValue; }
   };

   template <typename T>
   class unique_ptr<T[]> {
      T *fValue;
   public:
     VECPHYS_CUDA_HEADER_BOTH
     unique_ptr(T *in) : fValue(in) {}

     VECPHYS_CUDA_HEADER_BOTH
     ~unique_ptr() { delete [] fValue; }

     VECPHYS_CUDA_HEADER_BOTH
     T &operator[](size_t idx) { return fValue[idx]; }
   };
#endif

VECPHYS_GLOBAL int kAlignmentBoundary = 32;
VECPHYS_GLOBAL Precision kPi = 3.14159265358979323846;
VECPHYS_GLOBAL Precision kTwoPi = 2.*kPi;
VECPHYS_GLOBAL Precision kTwoPiInv = 1./kTwoPi;
VECPHYS_GLOBAL Precision kDegToRad = kPi/180.;
VECPHYS_GLOBAL Precision kRadToDeg = 180./kPi;
VECPHYS_GLOBAL Precision kInfinity =
#ifndef VECPHYS_NVCC
    std::numeric_limits<Precision>::infinity();
#else
    INFINITY;
#endif
VECPHYS_GLOBAL Precision kEpsilon =
#ifndef VECPHYS_NVCC
    std::numeric_limits<Precision>::epsilon();
#elif VECPHYS_FLOAT_PRECISION
    FLT_EPSILON;
#else
    DBL_EPSILON;
#endif
VECPHYS_GLOBAL Precision kMinimum =
#ifndef VECPHYS_NVCC
    std::numeric_limits<Precision>::min();
#elif VECPHYS_FLOAT_PRECISION
    FLT_MIN;
#else
    DBL_MIN;
#endif
VECPHYS_GLOBAL Precision kMaximum =
#ifndef VECPHYS_NVCC
    std::numeric_limits<Precision>::max();
#elif VECPHYS_FLOAT_PRECISION
    FLT_MAX;
#else
    DBL_MAX;
#endif
VECPHYS_GLOBAL Precision kTiny = 1e-30;
VECPHYS_GLOBAL Precision kTolerance = 1e-12;
VECPHYS_GLOBAL Precision kRadTolerance = 1e-12;
VECPHYS_GLOBAL Precision kAngTolerance = 1e-12;

VECPHYS_GLOBAL Precision kHalfTolerance = 0.5*kTolerance;
VECPHYS_GLOBAL Precision kToleranceSquared = kTolerance*kTolerance;

namespace EInside {
VECPHYS_GLOBAL vecphys::Inside_t kInside = 0;
VECPHYS_GLOBAL vecphys::Inside_t kSurface = 1;
VECPHYS_GLOBAL vecphys::Inside_t kOutside = 2;
}

typedef int RotationCode;
typedef int TranslationCode;
namespace rotation {
enum RotationId { kGeneric = -1, kDiagonal = 0x111, kIdentity = 0x200 };
}
namespace translation {
enum TranslationId { kGeneric = -1, kIdentity = 0 };
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
void Assert(const bool condition, char const *const message) {
#ifndef VECPHYS_NVCC
  if (!condition) {
    printf("Assertion failed: %s", message);
    abort();
  }
#else
  if (!condition) printf("Assertion failed: %s", message);
#endif
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
void Assert(const bool condition) {
  Assert(condition, "");
}

namespace details {
   template <typename DataType, typename Target> struct UseIfSameType { 
      VECPHYS_CUDA_HEADER_BOTH
      static Target const *Get(DataType*) { return nullptr; }
   };
   template <typename DataType> struct UseIfSameType<DataType,DataType> {
      VECPHYS_CUDA_HEADER_BOTH
      static DataType const *Get(DataType *ptr) { return ptr; }
   };
}

} } // End global namespace

#endif // VECPHYS_BASE_GLOBAL_H_
