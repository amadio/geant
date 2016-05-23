#ifndef VECPHYS_RANDOM_H
#define VECPHYS_RANDOM_H

#include <VecCore/VecCore>

#if !defined(VECCORE_NVCC)
#include <random>
#include <type_traits>
#else
#include <cuda.h>
#include <curand_kernel.h>
#endif

namespace vecCore {

#ifdef VECCORE_NVCC
namespace cuda {
static __shared__ curandState CudaRNGState;
}
#endif

template <typename T>
VECCORE_FORCE_INLINE VECCORE_CUDA_HOST_DEVICE T UniformRandom(void * /*States*/, int * /*ThreadId*/) {
#ifndef VECCORE_NVCC_DEVICE
  unsigned short xsubi[3];
  return erand48(xsubi);
#else
#if 1
  return curand_uniform(&cuda::CudaRNGState);
#else
  // for some reason, this implementation doesn't seem to work
  curandState *state_ptr = (curandState *)States;
  curandState localState = state_ptr[*ThreadId];
  T rng = curand_uniform(&localState);
  state_ptr[*ThreadId] = localState;
  return rng;
#endif
#endif
}

#if !defined(VECCORE_NVCC) && defined(VECCORE_ENABLE_VC)
#include <Vc/Vc>

template <> VECCORE_FORCE_INLINE Vc::Vector<float> UniformRandom<Vc::Vector<float>>(void *, int *) {
  return Vc::Vector<float>::Random();
}

template <> VECCORE_FORCE_INLINE Vc::Vector<double> UniformRandom<Vc::Vector<double>>(void *, int *) {
  return Vc::Vector<double>::Random();
}
#endif
}

#endif
