/// \file vc/backend.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)
//
// temporary clone for GUTracking based on VecGeom - syjun
//
#ifndef VECPHYS_SCALARBACKEND_H
#define VECPHYS_SCALARBACKEND_H

#include "base/Global.h"

#include <algorithm>
#include <cstring>
#include <stdlib.h>

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

struct kScalar {

  typedef int             int;
  typedef Precision       double;
  typedef bool            Bool_t;
  typedef int             Index_t; // the type of indices

  const static bool kTrue = true;
  const static bool kFalse = false;
  const static bool early_returns = true;

  const static int kSize = 1;
#ifdef VECPHYS_STD_CXX11
  constexpr static Precision kOne = 1.0;
  constexpr static Precision kZero = 0.0;
#else
  const static Precision kOne = 1.0;
  const static Precision kZero = 0.0;
#endif

  template <class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  static VECPHYS_CONSTEXPR_RETURN bool IsEqual() { return false; }

  VECPHYS_CUDA_HEADER_BOTH
  VECPHYS_INLINE
  static Precision Convert(Precision const &input) { return input; }
};

template <>
VECPHYS_CUDA_HEADER_BOTH
inline VECPHYS_CONSTEXPR_RETURN bool kScalar::IsEqual<kScalar>() {
  return true;
}

typedef kScalar::int    ScalarInt;
typedef kScalar::double ScalarDouble;
typedef kScalar::Bool_t   ScalarBool;

template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
void CondAssign(const bool cond,
                Type const &thenval, Type const &elseval, Type *const output) {
  *output = (cond) ? thenval : elseval;
}

template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
void MaskedAssign(const bool cond,
                  Type const &thenval, Type *const output) {
  *output = (cond) ? thenval : *output;
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
bool IsFull(bool const &cond){
    return cond;
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
bool Any(bool const &cond) {
  return cond;
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
bool IsEmpty(bool const &cond){
    return !cond;
}

template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Type Pow(Type const &x, Type arg) {
   return pow(x,arg);
}

template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Type Pow(Type const &x, int arg) {
   return pow(x,arg);
}

template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Type Abs(const Type val) {
  return fabs(val);
}

template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Type Sqrt(const Type val) {
  return std::sqrt(val);
}

template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Type Pow(const Type val1, const Type val2) {
  return std::pow(val1, val2);
}

template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Type Cbrt(const Type val1) {
  return cbrt(val1);
}


template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Type ATan2(const Type y, const Type x) {
  if (x != 0) return  std::atan2(y, x);
  if (y >  0) return  kPi / 2;
  if (y <  0) return -kPi / 2;
  return  0;
}

template <typename T>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
T Min(T const &val1, T const &val2) {
#ifndef VECPHYS_NVCC_DEVICE
  return std::min(val1, val2);
#else
  return val1 < val2 ? val1 : val2;
#endif
}

template <typename T>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
T Max(T const &val1, T const &val2) {
#ifndef VECPHYS_NVCC_DEVICE
  return std::max(val1, val2);
#else
  return val1 > val2 ? val1 : val2;
#endif
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Precision sin(const Precision radians) {
  return std::sin(radians);
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Precision cos(const Precision radians) {
  return std::cos(radians);
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Precision tan(const Precision radians) {
  return std::tan(radians);
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Precision Log(const Precision x) {
  return std::log(x);
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Precision Exp(const Precision x) {
  return std::exp(x);
}

VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Precision Floor( Precision val ){
    return std::floor(val);
}

template <typename Type>
VECPHYS_CUDA_HEADER_BOTH
VECPHYS_INLINE
Precision UniformRandom(Random_t* states, int tid){
#if __CUDA_ARCH__
  curandState localState = states[tid];
  Precision ran = curand_uniform(&localState);
  states[tid] = localState;
  return ran;
#else
  return (double)rand()/RAND_MAX;
#endif
}

} } // End global namespace

#endif // VECPHYS_SCALARBACKEND_H
