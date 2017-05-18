#ifndef VECRNG_H
#define VECRNG_H 1

/**
 * VecRNG: The base class of SIMD/SIMT random number generators
 *
 * Requirements :
 * 1) DerivedT  : A pseudo-random number generator with multiple streams
 * 2) BackendT  : Scalar, Vector, Cuda
 * 3) RandomT   : A templated struct of DerivedT states with BackendT 
 */

#include "base/VecPhys.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

template <typename DerivedT, typename BackendT, typename RandomT>
class VecRNG {

protected:
  // Use *this to access data members in the derived class
  RandomT *fState;

public:

  VECCORE_ATT_HOST_DEVICE
  VecRNG() { fState = new RandomT; }

  VECCORE_ATT_HOST_DEVICE
  ~VecRNG() { delete fState; }

  VECCORE_ATT_HOST_DEVICE
  VecRNG(const VecRNG &rng) = default;

  // Static interfaces (Required methods)
  
  // Initialization for SIMD
  VECCORE_ATT_HOST 
  void Initialize() 
  { static_cast<DerivedT *>(this)->template Initialize<BackendT>(); }

  // Initialization for SIMT
  VECCORE_ATT_HOST 
  void Initialize(RandomT *states, int blocks, int threads)
  { static_cast<DerivedT *>(this)->template Initialize<BackendT>(states,blocks,threads); }

  // Return Backend::Double_v of random numbers in [0,1)
  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE 
  typename Backend::Double_v Uniform() 
  { 
    return static_cast<DerivedT *>(this)->template Kernel<Backend>(this->fState); 
  }

  // Generate random numbers based on a given state
  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE 
  typename Backend::Double_v Uniform(RandomT *state) 
  { 
    return static_cast<DerivedT *>(this)->template Kernel<Backend>(state); 
  }

  // Auxiliary methods 

  VECCORE_ATT_HOST_DEVICE 
  void SetState(RandomT *state) { fState = state; }

  VECCORE_ATT_HOST_DEVICE 
  RandomT* GetState() const { return fState; }

  VECCORE_ATT_HOST void PrintState() const { static_cast<DerivedT *>(this)->PrintState(); }

  //Common methods

  // Returns an array of random numbers of the type Backend::Double_v
  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE 
  void Array(const size_t nsize, typename Backend::Double_v *array);

  // Flat distribution in [min,max)
  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v Flat(typename Backend::Double_v min, 
                                                           typename Backend::Double_v max) 
  {
    return min+(max-min)*static_cast<DerivedT *>(this)-> template Uniform<Backend>();
  }

  // Flat distribution in [min,max] with a state
  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v Flat(RandomT *state,
                                                          typename Backend::Double_v min, 
                                                          typename Backend::Double_v max) 
  {
    return min+(max-min)*static_cast<DerivedT *>(this)-> template Uniform<Backend>(state);
  }

  // Exponential deviates: exp(-x/tau)
  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v Exp(typename Backend::Double_v tau);

  // Exponential deviates with a state
  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v Exp(RandomT *state,
                                                         typename Backend::Double_v tau);

  // Gaussin deviates: 1/(2*pi*sigma^2)*exp[-(x-mean)^2/sigma^2]
  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v Gauss(typename Backend::Double_v mean, 
                                                           typename Backend::Double_v sigma);

  // Gaussin deviates with a state
  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v Gauss(RandomT *state,
                                                           typename Backend::Double_v mean, 
                                                           typename Backend::Double_v sigma);

  // @syj add functions to generate additional random distributions
  // (Binomial, Poisson, Landau and etc)

};

// Implementation

// Common Methods

// Returns an array of random numbers of Backend::Double_v
template <typename DerivedT, typename BackendT, typename RandomT>
template <typename Backend>
VECCORE_ATT_HOST_DEVICE void
VecRNG<DerivedT, BackendT, RandomT>::Array(const size_t nsize, typename Backend::Double_v *array)
{
  using Double_v = typename Backend::Double_v;
  for (size_t i = 0; i < nsize ; ++i) {
    Double_v u01 = static_cast<DerivedT *>(this)-> template Uniform<Backend>();    
    array[i] = u01;
  }
}
 
// Exponential deviates: exp(-x/tau)
template <typename DerivedT, typename BackendT, typename RandomT>
template <typename Backend>
VECCORE_ATT_HOST_DEVICE typename Backend::Double_v
VecRNG<DerivedT, BackendT, RandomT>::Exp(typename Backend::Double_v tau)
{
  using Double_v = typename Backend::Double_v;

  Double_v u01 = static_cast<DerivedT *>(this)-> template Uniform<Backend>();
  //@syjun check for zero 
  return -tau*math::Log(u01);
}

// Exponential deviates with state
template <typename DerivedT, typename BackendT, typename RandomT>
template <typename Backend>
VECCORE_ATT_HOST_DEVICE typename Backend::Double_v
VecRNG<DerivedT, BackendT, RandomT>::Exp(RandomT *state, typename Backend::Double_v tau)
{
  // Exp with state
  using Double_v = typename Backend::Double_v;

  Double_v u01 = static_cast<DerivedT *>(this)-> template Uniform<Backend>(state);
  return -tau*math::Log(u01);
}

// Gaussian deviates with state
template <typename DerivedT, typename BackendT, typename RandomT>
template <typename Backend>
VECCORE_ATT_HOST_DEVICE typename Backend::Double_v
VecRNG<DerivedT, BackendT, RandomT>::Gauss(RandomT *state, typename Backend::Double_v mean, 
                                           typename Backend::Double_v sigma)
{
  // Gauss with state
  using Double_v = typename Backend::Double_v;

  Double_v u1= static_cast<DerivedT *>(this)-> template Uniform<Backend>(state);
  Double_v u2= static_cast<DerivedT *>(this)-> template Uniform<Backend>(state)*(2.0*M_PI);
  Double_v normal =  math::Sqrt(-2.0*math::Log(u1))*math::Cos(u2);
  return mean+sigma*normal;
}

// Gaussian deviates with (mean, sigma):  1/(2*pi*sigma^2)*exp[-(x-mean)^2/sigma^2]
template <typename DerivedT, typename BackendT, typename RandomT>
template <typename Backend>
VECCORE_ATT_HOST_DEVICE typename Backend::Double_v
VecRNG<DerivedT, BackendT, RandomT>::Gauss(typename Backend::Double_v mean, typename Backend::Double_v sigma)
{
  // Using Box/Muller - use just one
  // normal1 = sqrt(-2*log(u1))*cos(2*pi*u2)
  // normal2 = sqrt(-2*log(u1))*sin(2*pi*u2)

  using Double_v = typename Backend::Double_v;

  Double_v u1= static_cast<DerivedT *>(this)-> template Uniform<Backend>();
  Double_v u2= static_cast<DerivedT *>(this)-> template Uniform<Backend>()*(2.0*M_PI);
  Double_v normal =  math::Sqrt(-2.0*math::Log(u1))*math::Cos(u2);
  return mean+sigma*normal;
}
 
} // end namespace impl
} // end namespace vecphys

#endif
