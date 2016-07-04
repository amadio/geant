#ifndef VECPHYS_BASE_GLOBAL_H_
#define VECPHYS_BASE_GLOBAL_H_

#define VECPHYS

#include <cstdio>

#include <VecCore/VecCore>

namespace vecphys {
using namespace vecCore;

using ScalarBackend = backend::Scalar;

#ifdef VECCORE_NVCC
using VectorBackend = backend::Scalar;
#else
using VectorBackend = backend::VcVector;
#endif

// scalar types

using Size_t = vecCore::Size_s;

using Real_t = vecCore::Real_s;
using Float_t = vecCore::Float_s;
using Double_t = vecCore::Double_s;

using Int_t = vecCore::Int_s;
using Int16_t = vecCore::Int16_s;
using Int32_t = vecCore::Int32_s;
using Int64_t = vecCore::Int64_s;

using UInt_t = vecCore::UInt_s;
using UInt16_t = vecCore::UInt16_s;
using UInt32_t = vecCore::UInt32_s;
using UInt64_t = vecCore::UInt64_s;

// vector types

using Real_v = typename VectorBackend::Real_v;
using Float_v = typename VectorBackend::Float_v;
using Double_v = typename VectorBackend::Double_v;

using Int_v = typename VectorBackend::Int_v;
using Int16_v = typename VectorBackend::Int16_v;
using Int32_v = typename VectorBackend::Int32_v;
using Int64_v = typename VectorBackend::Int64_v;

using UInt_v = typename VectorBackend::UInt_v;
using UInt16_v = typename VectorBackend::UInt16_v;
using UInt32_v = typename VectorBackend::UInt32_v;
using UInt64_v = typename VectorBackend::UInt64_v;
}

#ifdef VECCORE_NVCC
#include <cuda.h>
#include <curand_kernel.h>
typedef curandState Random_t;
#else
typedef int Random_t;
#endif

#include <Random.h>

#define VECPHYS_NAMESPACE ::vecphys

#if defined(VECCORE_NVCC)
#define VECPHYS_IMPL_NAMESPACE cuda
#else
#define VECPHYS_IMPL_NAMESPACE cxx
#endif

#ifdef VECCORE_NVCC_DEVICE
#define VECPHYS_GLOBAL static __constant__ const
#else
#define VECPHYS_GLOBAL constexpr
#endif

#define CLHEP vecphys // use CLHEP SystemOfUnits.h and PhysicalConstants.h

#if defined(VECCORE_NVCC)
#define VECPHYS_HOST_FORWARD_DECLARE(X)                                                                                \
  namespace cxx {                                                                                                      \
  X                                                                                                                    \
  }
#define VECPHYS_DEVICE_FORWARD_DECLARE(X)
#else
// Not compiling with NVCC
namespace vecphys {
template <typename DataType>
struct kCudaType;
template <typename DataType>
using CudaType_t = typename kCudaType<DataType>::type_t;
template <>
struct kCudaType<float> {
  using type_t = float;
};
template <>
struct kCudaType<double> {
  using type_t = double;
};
template <>
struct kCudaType<int> {
  using type_t = int;
};
}

#define VECPHYS_HOST_FORWARD_DECLARE(X)
#define VECPHYS_DEVICE_FORWARD_DECLARE(X)                                                                              \
  namespace cuda {                                                                                                     \
  X                                                                                                                    \
  }
#endif

#endif // VECPHYS_BASE_GLOBAL_H_
