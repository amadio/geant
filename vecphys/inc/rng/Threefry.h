#ifndef Threefry_H
#define Threefry_H 1

/**
  This is a SIMD/SIMT implementation of Threefry (4x32-20) 
  based on threefry.h of Random123 (version 1.09)
  
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

#include "R123.h"
#include <limits.h>

// struct Threefry_t (random state of Threefry-4x32-20)

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

template <typename BackendT>
struct Threefry_t {
  unsigned int index;
  R123::array4_t<BackendT> ctr;
  R123::array4_t<BackendT> key;
  R123::array4_t<BackendT> ukey;
}; 

//class Threefry<BackendT>

template <typename BackendT>
class Threefry : public VecRNG<Threefry<BackendT>, BackendT, Threefry_t<BackendT> > {

private:
  static long long fSeed;

public:

  VECCORE_ATT_HOST_DEVICE
  Threefry() {} 

  VECCORE_ATT_HOST_DEVICE
  ~Threefry() {}

  VECCORE_ATT_HOST_DEVICE
  Threefry(const Threefry &rng); 

  // Static methods

  inline VECCORE_ATT_HOST void Initialize();

  //Initialize a set of states of which size is equivalent to blocks*threads
  inline VECCORE_ATT_HOST void Initialize(Threefry_t<BackendT> *states, int blocks, int threads);

  // Returns pRNG<BackendT> between 0 and 1 (excluding the end points).
  template <typename Backend>
  inline VECCORE_ATT_HOST_DEVICE typename Backend::Double_v Kernel(Threefry_t<BackendT>& state);

  // Auxiliary methods

  VECCORE_ATT_HOST_DEVICE void SetSeed(long long seed) { fSeed = seed; }

  VECCORE_ATT_HOST_DEVICE long long GetSeed() const { return fSeed; }

  VECCORE_ATT_HOST void PrintState() const;

private:
  // the mother is friend of this
  friend class VecRNG<Threefry<BackendT>, BackendT, Threefry_t<BackendT> >;

  // Set the stream to the next stream/substream.
  VECCORE_ATT_HOST 
  void SetNextStream ();

  VECCORE_ATT_HOST 
  void SetNextSubstream ();

  // Increase counter
  VECCORE_ATT_HOST_DEVICE 
  void IncreaseCounter(Threefry_t<BackendT> *state);

  // Threefry utility methods
  VECCORE_ATT_HOST_DEVICE
  inline  typename BackendT::UInt32_v RotL_32(typename BackendT::UInt32_v x, unsigned int N);

  VECCORE_ATT_HOST_DEVICE
  inline  typename BackendT::UInt64_v RotL_64(typename BackendT::UInt64_v x, unsigned int N);

  VECCORE_ATT_HOST_DEVICE
  void BijectAndShuffle(R123::array4_t<BackendT> X, R123::array_t<BackendT,5> ks, 
                        unsigned int start, unsigned int index);

  VECCORE_ATT_HOST_DEVICE 
  void Gen(R123::array4_t<BackendT> in, R123::array4_t<BackendT> k, R123::array4_t<BackendT> X);

};

// The default seed of Threefry
template <class BackendT> long long Threefry<BackendT>::fSeed = 12345;     

//
// Class Implementation
//  

// Copy constructor
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE
Threefry<BackendT>::Threefry(const Threefry<BackendT> &rng) 
  : VecRNG<Threefry<BackendT>, BackendT, Threefry_t<BackendT> >()
{
  this->fState->index = rng.fState->index;

  for(size_t i = 0 ; i < 4 ; ++i) {
    this->fState->ctr[i] = rng.fState->ctr[i];
    this->fState->key[i] = rng.fState->key[i];
    this->fState->ukey[i]= rng.fState->ukey[i];
  }
  fSeed = rng.fSeed; 
}

// Set a new set of keys for the vector of next streams using the unique seed
template <typename BackendT>
VECCORE_ATT_HOST void Threefry<BackendT>::SetNextStream ()
{
  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) { 
    this->fState->key[0][i] = (unsigned int)(fSeed);
    this->fState->key[1][i] = (unsigned int)(fSeed>>32);
    ++fSeed; 
  }
}

// Scalar Specialization of SetNextStream
template <>
VECCORE_ATT_HOST void Threefry<ScalarBackend>::SetNextStream ()
{
  this->fState->key[0] = (unsigned int)(fSeed);
  this->fState->key[1] = (unsigned int)(fSeed>>32);
  ++fSeed; 
}

// Reset the current stream to the next substream - skipahead 
template <typename BackendT>
VECCORE_ATT_HOST void Threefry<BackendT>::SetNextSubstream ()
{  
  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) { 
    unsigned int nlo = (unsigned int)(R123::THREEFRY_SKIP_AHEAD);
    unsigned int nhi = (unsigned int)(R123::THREEFRY_SKIP_AHEAD>>32);

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
VECCORE_ATT_HOST void Threefry<ScalarBackend>::SetNextSubstream ()
{
  unsigned int nlo = (unsigned int)(R123::THREEFRY_SKIP_AHEAD);
  unsigned int nhi = (unsigned int)(R123::THREEFRY_SKIP_AHEAD>>32);

  this->fState->ctr[0] += nlo;
  if( this->fState->ctr[0] < nlo ) nhi++;
  this->fState->ctr[1] += nhi;
  if(nhi <= this->fState->ctr[1]) return;
  if(++this->fState->ctr[2]) return;
  ++this->fState->ctr[3];
}

//
// set the seed for each stream of SIMD to the starting state of each substream
//
template <class BackendT>
VECCORE_ATT_HOST void Threefry<BackendT>::Initialize()
{
  //set initial counter and key
  this->fState->index = 0;
  SetNextStream();
}

// Specialization of Initialize for SIMT
template <>
VECCORE_ATT_HOST void Threefry<ScalarBackend>::Initialize(Threefry_t<ScalarBackend> *states, 
                                                          int blocks, int threads)
{
  
  Threefry_t<ScalarBackend>* hstates 
    = (Threefry_t<ScalarBackend> *) malloc (blocks*threads*sizeof(Threefry_t<ScalarBackend>));

  unsigned int nthreads = blocks * threads;
  for (unsigned int tid = 0 ; tid < nthreads ; ++tid) {
    //initialize initial seed/state by the unique tid number
    hstates[tid].index = 0; 
    fSeed = tid;
    SetNextStream();
    //assume that elments of ctr and ukey are set to zero
  }

#ifdef VECCORE_CUDA
  cudaMemcpy(states, hstates, nthreads*sizeof(Threefry_t<ScalarBackend>), cudaMemcpyHostToDevice);
#else
  memcpy(states, hstates, nthreads*sizeof(Threefry_t<ScalarBackend>));
#endif
  free(hstates);
}

// Increase counter of each element of the counter (ctr) vector
template <typename BackendT>
VECCORE_ATT_HOST void Threefry<BackendT>::IncreaseCounter(Threefry_t<BackendT> *state)
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
VECCORE_ATT_HOST void Threefry<ScalarBackend>::IncreaseCounter(Threefry_t<ScalarBackend> *state)
{
  if( ++state->ctr[0]) return;
  if( ++state->ctr[1]) return;
  if( ++state->ctr[2]) return;
  ++state->ctr[3];
}

// Print information of the current state
template <typename BackendT>
VECCORE_ATT_HOST void Threefry<BackendT>::PrintState() const
{
  std::cout << "index = " << this->fState->index << std::endl;
  for(size_t i = 0 ; i < 4 ; ++i) {
    std::cout << "key[" << i << "] = " <<  this->fState->key[i] << std::endl;
  }
}

// Kernel to generate a vector(scalar) of next random number(s) 
template <class BackendT>
template <class Backend>
VECCORE_ATT_HOST_DEVICE typename Backend::Double_v Threefry<BackendT>::Kernel(Threefry_t<BackendT>& state)
{

  /*
  if(last_elem == 0){
    // jam n into the high bits of c
    const size_t W = std::numeric_limits<result_type>::digits;
    ctr_type c = c0;
    c[c0.size()-1] |= n<<(W-BITS);
    rdata = b(c,k);
    n++;
    last_elem = rdata.size();
  }
  return rdata[--last_elem];
  */
  using Double_v = typename Backend::Double_v;
  Double_v u(0.0);

  if(state.index == 0 ) {
    //get a new state and generate 128 bits of pseudo randomness 
    Gen(state.ctr,state.key,state.ukey);

    //construct 8xUInt32 (ukey) to 4xUInt64, then convert to Double_v using UINT64_MAX
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

// Threefry utility methods
template <class BackendT>
inline VECCORE_ATT_HOST_DEVICE typename BackendT::UInt32_v 
Threefry<BackendT>::RotL_32(typename BackendT::UInt32_v x, unsigned int N)
{
  return (x << (N & 31)) | (x >> ((32-N) & 31));
}

template <class BackendT>
inline VECCORE_ATT_HOST_DEVICE typename BackendT::UInt64_v 
Threefry<BackendT>::RotL_64(typename BackendT::UInt64_v x, unsigned int N)
{
  return (x << (N & 63)) | (x >> ((64-N) & 63));
}

template <class BackendT>
VECCORE_ATT_HOST_DEVICE void 
Threefry<BackendT>::BijectAndShuffle(R123::array4_t<BackendT> X, R123::array_t<BackendT,5> key, 
                                     unsigned int start, unsigned int index) 
{
  X[0] += X[1]; X[1] = RotL_32(X[1],R123::R_32x4[start+0][0]); X[1] ^= X[0];
  X[2] += X[3]; X[3] = RotL_32(X[3],R123::R_32x4[start+0][1]); X[3] ^= X[2];

  X[0] += X[3]; X[3] = RotL_32(X[3],R123::R_32x4[start+1][0]); X[3] ^= X[0];
  X[2] += X[1]; X[1] = RotL_32(X[1],R123::R_32x4[start+1][1]); X[1] ^= X[2];

  X[0] += X[1]; X[1] = RotL_32(X[1],R123::R_32x4[start+2][0]); X[1] ^= X[0];
  X[2] += X[3]; X[3] = RotL_32(X[3],R123::R_32x4[start+2][1]); X[3] ^= X[2];

  X[0] += X[3]; X[3] = RotL_32(X[3],R123::R_32x4[start+3][0]); X[3] ^= X[0];
  X[2] += X[1]; X[1] = RotL_32(X[1],R123::R_32x4[start+3][1]); X[1] ^= X[2];

  // InjectKey (r=index) 			
  for(size_t ir = 0 ; ir < 4 ; ++ir) {
    int kpos = (ir+index) % 5;
    X[ir] += key[kpos];
  }
  X[3] += index;     // X[WCNT4-1] += r  
}  

template <class BackendT>
VECCORE_ATT_HOST_DEVICE void Threefry<BackendT>::Gen(R123::array4_t<BackendT> in, R123::array4_t<BackendT> k, 
                                                     R123::array4_t<BackendT> X)
{
  using UInt32_v = typename BackendT::UInt32_v;                              

  UInt32_v ks[4+1];				
  ks[4] = 0x1BD11BDA;

  //  using UInt64_v = typename BackendT::UInt64_v;                              
  //  UInt64_v ks64[5];				
  //  ks64[4] = (0x1BD11BDA) + (((UInt64_v) (0xA9FC1A22)) << 32) ;

  //  UInt32_v ks32[5];				
  //  ks32[4] = (0x1BD11BDA) + (((UInt32_v) (0xA9FC1A22)) << 16) ;
  
  for (size_t i = 0 ; i < 4 ; ++i) {                                       
    ks[i] = k[i];
    X[i]  = in[i];
    ks[4] ^= k[i];
  }      

  // Insert initial key before round 0     
  X[0] += ks[0];
  X[1] += ks[1];
  X[2] += ks[2];
  X[3] += ks[3];				

  //round 1-4
  BijectAndShuffle(X,ks,0,1);

  //round 5-8
  BijectAndShuffle(X,ks,4,2);

  //round 9-12
  BijectAndShuffle(X,ks,0,3);

  //round 13-16
  BijectAndShuffle(X,ks,4,4);

  //round 17-20
  BijectAndShuffle(X,ks,0,5);
}

} // end namespace impl
} // end namespace vecphys

#endif
