#ifndef VECPHYS_RANDOM_H
#define VECPHYS_RANDOM_H

#include "base/VecPhys.h"
#include <VecCore/VecCore>

#if !defined(VECCORE_CUDA)
#include <random>
#include <type_traits>
#else
#include <cuda.h>
#include <curand_kernel.h>
#endif

namespace vecCore {

template <typename T>
VECCORE_FORCE_INLINE VECCORE_ATT_HOST_DEVICE T UniformRandom(Random_t *states, int threadId)
{
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  (void)states;
  (void)threadId;
  return ((double)rand()) / ((double)(RAND_MAX));
#else
  curandState localState = states[threadId];
  T rng = curand_uniform(&localState);
  states[threadId] = localState;
  return rng;
#endif
}

#if !defined(VECCORE_CUDA) && defined(VECCORE_ENABLE_VC)
#include <Vc/Vc>

template <>
VECCORE_FORCE_INLINE Vc::Vector<float> UniformRandom<Vc::Vector<float>>(Random_t *, int)
{
  return Vc::Vector<float>::Random();
}

template <>
VECCORE_FORCE_INLINE Vc::Vector<double> UniformRandom<Vc::Vector<double>>(Random_t *, int)
{
  return Vc::Vector<double>::Random();
}
#endif
}

#endif
