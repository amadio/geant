#ifndef VECRNG_H
#define VECRNG_H 1

/**
 * VecRNG: The base class of pRNGs portable for both SIMD and SIMT
 *
 * Requirements:
 * 1) pRNG     : A pseudo-random number generator with multiple streams
 * 2) Type     : Scalar, Vector, Cuda
 * 3) Random_t : A templated struct of pRNG states with Type 
 */

#include "base/VecPhys.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

template <typename pRNG, typename Type, typename Random_t>
class VecRNG {

protected:
  // Use *this to access data members in the derived class
  Random_t *fState;

public:

  VECCORE_ATT_HOST_DEVICE
  VecRNG() { fState = new Random_t; }

  VECCORE_ATT_HOST_DEVICE
  ~VecRNG() { delete fState; }

  VECCORE_ATT_HOST_DEVICE
  VecRNG(const VecRNG &rng) = default;

  // Static interfaces (Required methods)
  
  // Initialization for SIMD
  VECCORE_ATT_HOST 
  void Initialize() 
  { static_cast<pRNG *>(this)->template Initialize<Type>(); }

  // Initialization for SIMT
  VECCORE_ATT_HOST 
  void Initialize(Random_t *states, int blocks, int threads)
  { static_cast<pRNG *>(this)->template Initialize<Type>(states,blocks,threads); }

  // Return RetType::Double_v of random numbers in [0,1)
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE 
  typename RetType::Double_v Uniform() { return static_cast<pRNG *>(this)->template Kernel<RetType>(this->fState); }

  // Generate random numbers based on a given state
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE 
  typename RetType::Double_v Uniform(Random_t *state) 
  { 
    return static_cast<pRNG *>(this)->template Kernel<RetType>(state); 
  }

  // Auxiliary methods 

  VECCORE_ATT_HOST_DEVICE 
  void SetState(Random_t *state) { fState = state; }

  VECCORE_ATT_HOST_DEVICE 
  Random_t* GetState() const { return fState; }

  VECCORE_ATT_HOST void PrintState() const { static_cast<pRNG *>(this)->PrintState(); }

  //Common methods

  // Returns an array of random numbers of the type RetType::Double_v
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE 
  void Array(const size_t nsize, typename RetType::Double_v *array);

  // Flat distribution in [min,max] with a state
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE typename RetType::Double_v Flat(Random_t *state,
                                                          typename RetType::Double_v min, 
                                                          typename RetType::Double_v max);

  // Flat distribution in [min,max)
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE typename RetType::Double_v Flat(typename RetType::Double_v min, 
                                                          typename RetType::Double_v max);

  // Exponential deviates: exp(-x/tau)
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE typename RetType::Double_v Exp(typename RetType::Double_v tau);

  // Exponential deviates with a state
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE typename RetType::Double_v Exp(Random_t *state,
                                                         typename RetType::Double_v tau);

  // Gaussin deviates: 1/(2*pi*sigma^2)*exp[-(x-mean)^2/sigma^2]
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE typename RetType::Double_v Gauss(typename RetType::Double_v mean, 
                                                           typename RetType::Double_v sigma);

  // Gaussin deviates with a state
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE typename RetType::Double_v Gauss(Random_t *state,
                                                           typename RetType::Double_v mean, 
                                                           typename RetType::Double_v sigma);

  // @syj add functions to generate additional random distributions
  // (Binomial, Poisson, Landau and etc)

};

// Implementation

// Common Methods

// Returns an array of random numbers of RetType::Double_v
template <typename pRNG, typename Type, typename Random_t>
template <typename RetType>
VECCORE_ATT_HOST_DEVICE void
VecRNG<pRNG, Type, Random_t>::Array(const size_t nsize, typename RetType::Double_v *array)
{
  using Double_v = typename RetType::Double_v;
  for (size_t i = 0; i < nsize ; ++i) {
    Double_v u01 = static_cast<pRNG *>(this)-> template Uniform<RetType>();    
    array[i] = u01;
  }
}
 
// Flat distribution in [min, max] with a state
template <typename pRNG, typename Type, typename Random_t>
template <typename RetType>
VECCORE_ATT_HOST_DEVICE typename RetType::Double_v
VecRNG<pRNG, Type, Random_t>::Flat(Random_t *state, typename RetType::Double_v min, typename RetType::Double_v max)
{
  return min+(max-min)*static_cast<pRNG *>(this)-> template Uniform<RetType>(state);
}

// Flat distribution in [min, max] 
template <typename pRNG, typename Type, typename Random_t>
template <typename RetType>
VECCORE_ATT_HOST_DEVICE typename RetType::Double_v
VecRNG<pRNG, Type, Random_t>::Flat(typename RetType::Double_v min, typename RetType::Double_v max)
{
  return min+(max-min)*static_cast<pRNG *>(this)-> template Uniform<RetType>();
}

// Exponential deviates: exp(-x/tau)
template <typename pRNG, typename Type, typename Random_t>
template <typename RetType>
VECCORE_ATT_HOST_DEVICE typename RetType::Double_v
VecRNG<pRNG, Type, Random_t>::Exp(typename RetType::Double_v tau)
{
  using Double_v = typename RetType::Double_v;

  Double_v u01 = static_cast<pRNG *>(this)-> template Uniform<RetType>();
  //@syjun check for zero 
  return -tau*math::Log(u01);
}

// Exponential deviates with state
template <typename pRNG, typename Type, typename Random_t>
template <typename RetType>
VECCORE_ATT_HOST_DEVICE typename RetType::Double_v
VecRNG<pRNG, Type, Random_t>::Exp(Random_t *state, typename RetType::Double_v tau)
{
  // Exp with state
  using Double_v = typename RetType::Double_v;

  Double_v u01 = static_cast<pRNG *>(this)-> template Uniform<RetType>(state);
  return -tau*math::Log(u01);
}

// Gaussian deviates with state
template <typename pRNG, typename Type, typename Random_t>
template <typename RetType>
VECCORE_ATT_HOST_DEVICE typename RetType::Double_v
VecRNG<pRNG, Type, Random_t>::Gauss(Random_t *state, typename RetType::Double_v mean, typename RetType::Double_v sigma)
{
  // Gauss with state
  using Double_v = typename RetType::Double_v;

  Double_v u1= static_cast<pRNG *>(this)-> template Uniform<RetType>(state);
  Double_v u2= static_cast<pRNG *>(this)-> template Uniform<RetType>(state)*(2.0*M_PI);
  Double_v normal =  math::Sqrt(-2.0*math::Log(u1))*math::Cos(u2);
  return mean+sigma*normal;
}

// Gaussian deviates with (mean, sigma):  1/(2*pi*sigma^2)*exp[-(x-mean)^2/sigma^2]
template <typename pRNG, typename Type, typename Random_t>
template <typename RetType>
VECCORE_ATT_HOST_DEVICE typename RetType::Double_v
VecRNG<pRNG, Type, Random_t>::Gauss(typename RetType::Double_v mean, typename RetType::Double_v sigma)
{
  // Using Box/Muller - use just one
  // normal1 = sqrt(-2*log(u1))*cos(2*pi*u2)
  // normal2 = sqrt(-2*log(u1))*sin(2*pi*u2)

  using Double_v = typename RetType::Double_v;

  Double_v u1= static_cast<pRNG *>(this)-> template Uniform<RetType>();
  Double_v u2= static_cast<pRNG *>(this)-> template Uniform<RetType>()*(2.0*M_PI);
  Double_v normal =  math::Sqrt(-2.0*math::Log(u1))*math::Cos(u2);
  return mean+sigma*normal;
}
 
} // end namespace impl
} // end namespace vecphys

#endif
