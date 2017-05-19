#ifndef Philox_H
#define Philox_H 1

/**
  Philox: A SIMD/SIMT implementation of Philox (4x32-10) based on philox.h 
          of Random123 (version 1.09)
  
  Copyright 2010-2011, D. E. Shaw Research.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:
  
  * Redistributions of source code must retain the above copyright
    notice, this list of conditions, and the following disclaimer.
  
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions, and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  
  * Neither the name of D. E. Shaw Research nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "VecRNG.h"
#include "R123.h"

#include <limits.h>

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

// struct Philox (random state of Philox-4x32-10)

template <typename BackendT>
struct Philox_t {
  R123::array_t<BackendT,4> ctr;
  R123::array_t<BackendT,2> key;
  R123::array_t<BackendT,4> ukey;
  unsigned int index;
}; 

template <typename BackendT>
class Philox : public VecRNG<Philox<BackendT>, BackendT, Philox_t<BackendT> > {

private:
  static long long fSeed;

public:

  VECCORE_ATT_HOST_DEVICE
  Philox() {} 

  VECCORE_ATT_HOST_DEVICE
  ~Philox() {}

  VECCORE_ATT_HOST_DEVICE
  Philox(const Philox &rng); 

  // Static methods

  inline VECCORE_ATT_HOST void Initialize();

  //Initialize a set of states of which size is equivalent to blocks*threads
  inline VECCORE_ATT_HOST void Initialize(Philox_t<BackendT> *states, int blocks, int threads);

  // Returns pRNG<BackendT> between 0 and 1 (excluding the end points).
  template <typename Backend>
  inline VECCORE_ATT_HOST_DEVICE typename Backend::Double_v Kernel(Philox_t<BackendT>& state);

  // Auxiliary methods

  VECCORE_ATT_HOST_DEVICE void SetSeed(long long seed) { fSeed = seed; }

  VECCORE_ATT_HOST_DEVICE long long GetSeed() const { return fSeed; }

  VECCORE_ATT_HOST void PrintState() const;

private:
  // the mother is friend of this
  friend class VecRNG<Philox<BackendT> , BackendT, Philox_t<BackendT> >;


  // Set the stream to the next stream/substream.
  VECCORE_ATT_HOST void SetNextStream ();

  VECCORE_ATT_HOST void SetNextSubstream ();

  // Increase counter
  VECCORE_ATT_HOST_DEVICE void IncreaseCounter(Philox_t<BackendT> *state);
  
  // Philox utility methods
  VECCORE_ATT_HOST_DEVICE
  inline  typename BackendT::UInt32_v Mulhilo32(typename BackendT::UInt32_v a, typename BackendT::UInt32_v b, 
                                                typename BackendT::UInt32_v *hip);

  VECCORE_ATT_HOST_DEVICE inline void Philox4x32bumpkey(R123::array_t<BackendT,2> key);

  VECCORE_ATT_HOST_DEVICE void Philox4x32round(R123::array_t<BackendT,4> crt, 
					       R123::array_t<BackendT,2> key);
  VECCORE_ATT_HOST_DEVICE void Gen(R123::array_t<BackendT,4> ctr, R123::array_t<BackendT,2> key,
                                   R123::array_t<BackendT,4> out);

};

// The default seed of Philox
template <class BackendT> long long Philox<BackendT>::fSeed = 12345;

//
// Class Implementation
//  

// Copy constructor
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE
Philox<BackendT>::Philox(const Philox<BackendT> &rng) : VecRNG<Philox<BackendT>, BackendT, Philox_t<BackendT> >()
{
  this->fState->index = rng.fState->index;
  for(size_t i = 0 ; i < 4 ; ++i) {
    this->fState->ctr[i] = rng.fState->ctr[i];
    if(i < 2) this->fState->key[i] = rng.fState->key[i];
    this->fState->ukey[i]= rng.fState->ukey[i];
  }  
  fSeed = rng.fSeed; 
}

// Set a new set of keys for the vector of next streams using the unique seed
template <typename BackendT>
VECCORE_ATT_HOST void Philox<BackendT>::SetNextStream ()
{
  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) { 
    this->fState->key[0][i] = (unsigned int)(fSeed);
    this->fState->key[1][i] = (unsigned int)(fSeed>>32);
    ++fSeed; 
  }
}

// Scalar Specialization of SetNextStream
template <>
VECCORE_ATT_HOST void Philox<ScalarBackend>::SetNextStream ()
{
  this->fState->key[0] = (unsigned int)(fSeed);
  this->fState->key[1] = (unsigned int)(fSeed>>32);
  ++fSeed; 
}

// Reset the current stream to the next substream - skipahead 
template <typename BackendT>
VECCORE_ATT_HOST void Philox<BackendT>::SetNextSubstream ()
{
  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) { 
    unsigned int nlo = (unsigned int)(R123::PHILOX_SKIP_AHEAD);
    unsigned int nhi = (unsigned int)(R123::PHILOX_SKIP_AHEAD>>32);

    this->fState->ctr[0] += nlo;
    if( this->fState->ctr[0] < nlo ) nhi++;
    this->fState->ctr[1] += nhi;
    if(nhi <= this->fState->ctr[1]) continue;
    if(++this->fState->ctr[2]) continue;
    ++this->fState->ctr[3];
  }
}

// Scalar specialization of SetNextSubstream
template <>
VECCORE_ATT_HOST void Philox<ScalarBackend>::SetNextSubstream ()
{
  unsigned int nlo = (unsigned int)(R123::PHILOX_SKIP_AHEAD);
  unsigned int nhi = (unsigned int)(R123::PHILOX_SKIP_AHEAD>>32);

  this->fState->ctr[0] += nlo;
  if( this->fState->ctr[0] < nlo ) nhi++;
  this->fState->ctr[1] += nhi;
  if(nhi <= this->fState->ctr[1]) return;
  if(++this->fState->ctr[2]) return;
  ++this->fState->ctr[3];
}

// Set the seed for each stream of SIMD to the starting state of each substream
template <class BackendT>
inline VECCORE_ATT_HOST void Philox<BackendT>::Initialize()
{
  //set initial counter and key
  this->fState->index = 0;
  SetNextStream();
}

// Specialization of Initialize for SIMT
template <>
inline VECCORE_ATT_HOST void Philox<ScalarBackend>::Initialize(Philox_t<ScalarBackend> *states, 
                                                               int blocks, int threads)
{
  Philox_t<ScalarBackend>* hstates 
   = (Philox_t<ScalarBackend> *) malloc (blocks*threads*sizeof(Philox_t<ScalarBackend>));

  unsigned int nthreads = blocks * threads;

  for (unsigned int tid = 0 ; tid < nthreads ; ++tid) {
    //initialize initial seed/state by the unique tid number
    hstates[tid].index = 0; 
    fSeed = tid;
    SetNextStream();
    //assume that elments of ctr and ukey are set to zero
  }
#ifdef VECCORE_CUDA
  cudaMemcpy(states, hstates, nthreads*sizeof(Philox_t<ScalarBackend>), cudaMemcpyHostToDevice);
#else
  memcpy(states, hstates, nthreads*sizeof(Philox_t<ScalarBackend>));
#endif
  free(hstates);
}

// Increase counter of each element of the counter (ctr) vector
template <typename BackendT>
VECCORE_ATT_HOST void Philox<BackendT>::IncreaseCounter(Philox_t<BackendT> *state)
{
  size_t vsize = VectorSize<UInt32_v>();
  for(size_t iv = 0 ; iv < vsize ; ++iv) {
    if( ++state->ctr[0][iv]) continue;
    if( ++state->ctr[1][iv]) continue;
    if( ++state->ctr[2][iv]) continue;
    ++state->ctr[3][iv];
  }
}

// Increase counter of each element of the counter (ctr) vector
template <>
VECCORE_ATT_HOST void Philox<ScalarBackend>::IncreaseCounter(Philox_t<ScalarBackend> *state)
{
  if( ++state->ctr[0]) return;
  if( ++state->ctr[1]) return;
  if( ++state->ctr[2]) return;
  ++state->ctr[3];
}

// Print information of the current state
template <typename BackendT>
VECCORE_ATT_HOST void Philox<BackendT>::PrintState() const
{
  std::cout << "index = " << this->fState->index << std::endl;
  for(size_t i = 0 ; i < 2 ; ++i) {
    std::cout << "key[" << i << "] = " <<  this->fState->key[i] << std::endl;
  }
}

// Kernel to generate a vector(scalar) of next random number(s)
template <class BackendT>
template <class Backend>
inline VECCORE_ATT_HOST_DEVICE  typename Backend::Double_v Philox<BackendT>::Kernel(Philox_t<BackendT>& state)
{
  using Double_v = typename Backend::Double_v;
  Double_v u(0.0);

  if(state.index == 0 ) {
    //get a new state and generate 128 bits of pseudo randomness 
    Gen(state.ctr,state.key,state.ukey);

    //construct 8xUInt32 (ukey) to 4xUInt64, then convert to Double_v using UINT32_MAX
    u = static_cast<Double_v>( (state.ukey[state.index]) )/UINT32_MAX;
    
    //state index and increase counter
    ++(state.index);
    IncreaseCounter(&state);
  }
  else {  
    u = static_cast<Double_v>( (state.ukey[state.index]) )/UINT32_MAX; 
    ++state.index;
    if(state.index == 4) state.index = 0;
  }
  
  return u;
}

// Philox utility methods
template <class BackendT>
inline
VECCORE_ATT_HOST_DEVICE typename BackendT::UInt32_v Philox<BackendT>::Mulhilo32(typename BackendT::UInt32_v a,
									        typename BackendT::UInt32_v b, 
									        typename BackendT::UInt32_v *hip)
{
  using UInt32_v = typename BackendT::UInt32_v;                              
  using UInt64_v = typename BackendT::UInt64_v;                              

  //  UInt64_v product = ((UInt64_v)a)*((UInt64_v)b);        
  //  *hip = product >> 32;                                  
  //  return (UInt32_v)product; 

  //  @@@syj - we do not have Vc::UInt64_v for AVX/AVX2
  //  use sclar for now

  UInt64_v product;
  UInt32_v result;

  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) { 
    product[i] = ((uint64_t)a[i])*((uint64_t)b[i]);
    (*hip)[i] = product[i] >> 32;                                  
    result[i] = ((uint32_t)product[i]);
  }
  return result; 
}

template <>
inline
VECCORE_ATT_HOST_DEVICE ScalarBackend::UInt32_v Philox<ScalarBackend>::Mulhilo32(uint32_t a,
								                 uint32_t b, 
								                 uint32_t *hip)
{
#ifndef __CUDA_ARCH__
  uint64_t product = ((uint64_t)a)*((uint64_t)b);
  *hip = product >> 32;                                  
  return (uint32_t)product; 
#else
  *hip = __umulhi(a,b);
  return a*b;
#endif
}

template <class BackendT>
inline
VECCORE_ATT_HOST_DEVICE void Philox<BackendT>::Philox4x32bumpkey(R123::array_t<BackendT,2> key) {
  //  key[0] += (UInt32_v) R123::PHILOX_W32_0;                                       
  //  key[1] += (UInt32_v) R123::PHILOX_W32_1;                                       
  key[0] = key[0] + 0x9E3779B9;                                        
  key[1] = key[1] + 0xBB67AE85;                                       
}  

template <class BackendT>
VECCORE_ATT_HOST_DEVICE void Philox<BackendT>::Philox4x32round(R123::array_t<BackendT,4> ctr, 
                                                               R123::array_t<BackendT,2> key)
{
  using UInt32_v = typename BackendT::UInt32_v;                              

  UInt32_v hi0;                                                        
  UInt32_v hi1;                                                        

  //  UInt32_v lo0 = Mulhilo32(R123::PHILOX_M4x32_0, ctr[0], &hi0);            
  //  UInt32_v lo1 = Mulhilo32(R123::PHILOX_M4x32_1, ctr[2], &hi1);            
  UInt32_v lo0 = Mulhilo32( (UInt32_v)0xD2511F53, ctr[0], &hi0);            
  UInt32_v lo1 = Mulhilo32( (UInt32_v)0xCD9E8D57, ctr[2], &hi1);            

  //  UInt32_v out[4] = {hi1^ctr[1]^key[0], lo1, hi0^ctr[3]^key[1], lo0};    

  ctr[0] = hi1^ctr[1]^key[0];
  ctr[1] = lo1;
  ctr[2] = hi0^ctr[3]^key[1];
  ctr[3] = lo1;
}  

template <class BackendT>
VECCORE_ATT_HOST_DEVICE void Philox<BackendT>::Gen(R123::array_t<BackendT,4> ctr, R123::array_t<BackendT,2> key,
                                                   R123::array_t<BackendT,4> output)
{
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);   
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key);

  //  return ctr;                                                         
  for (int i=0;i < 4; i++) { 
    output[i] = ctr[i];                                                         
  }
}

} // end namespace impl
} // end namespace vecphys

#endif
