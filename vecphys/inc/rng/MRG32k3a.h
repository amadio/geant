#ifndef MRG32k3a_H
#define MRG32k3a_H 1

/**
 * MRG32k3a: A SIMD/SIMT implementation of MRG32k3a based on RngStream.h(cpp)
 *
 * RngStream is a class generating multiple streams of random numbers created
 * by Prof. Pierre L'Ecuyer, University of Montreal (lecuyer@iro.umontreal.ca) 
 * Original source codes of RngStream.h(cpp) is available at
 * http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c++/
 *
 * Relevant articles in which MRG32k3a and the package with multiple streams 
 * were proposed:
 *
 * P. L'Ecuyer, ``Good Parameter Sets for Combined Multiple Recursive Random
 * Number Generators'', Operations Research, 47, 1 (1999), 159--164.
 *
 * P. L'Ecuyer, R. Simard, E. J. Chen, and W. D. Kelton, ``An Objected-Oriented
 * Random-Number Package with Many Long Streams and Substreams'', Operations
 * Research, 50, 6 (2002), 1073--1075
 */

#include "VecRNG.h"
#include "MRG.h"

#include <iostream>

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

// struct MRG32k3a_t (random state of MRG32k3a)

template <typename Type>
struct MRG32k3a_t {
  typename Type::Double_v fCg[MRG::vsize];
}; 

//class MRG32k3a<Type>

template <typename Type>
class MRG32k3a : public VecRNG<MRG32k3a<Type>, Type, MRG32k3a_t<Type> > {

private:
  static Real_t fSeed[MRG::vsize];
  typename Type::Double_v Bg[MRG::vsize];

  // Information on a stream: The arrays {Cg, Bg, Ig} (from the RngStream) 
  // contain the current state of the stream, the starting state of the current
  // SubStream, and the starting state of the stream (not used in this class). 
  // The next seed will be the seed of the next declared RngStream. 

public:

  VECCORE_ATT_HOST_DEVICE
  MRG32k3a() {}

  VECCORE_ATT_HOST_DEVICE
  ~MRG32k3a() {}

  VECCORE_ATT_HOST_DEVICE
  MRG32k3a(const MRG32k3a &rng); 

  // Static methods

  VECCORE_ATT_HOST void Initialize() { SetNextStream(); }

  VECCORE_ATT_HOST void Initialize(MRG32k3a_t<Type> *states, int blocks, int threads);

  // Returns pRNG<Type> between 0 and 1 
  template <typename RetType>
  VECCORE_ATT_HOST_DEVICE typename RetType::Double_v Kernel(MRG32k3a_t<Type> *state);

  // Auxiliary methods

  VECCORE_ATT_HOST_DEVICE void SetSeed(Real_t seed[MRG::vsize]);

  VECCORE_ATT_HOST_DEVICE Real_t* GetSeed() const { return fSeed; }

  VECCORE_ATT_HOST void PrintState() const ;

private:
  
  // the mother is friend of this
  friend class VecRNG<MRG32k3a<Type>, Type, MRG32k3a_t<Type> >;

  // Set Stream to NextStream/NextSubStream.
  VECCORE_ATT_HOST void SetNextStream ();

  VECCORE_ATT_HOST void SetNextSubstream ();

  // MRG32k3a utility methods
  VECCORE_ATT_HOST double MultModM (double a, double s, double c, double m);

  VECCORE_ATT_HOST void MatVecModM (const double A[MRG::ndim][MRG::ndim],
				    const double s[MRG::ndim], double v[MRG::ndim], double m);
};

// The default seed of MRG32k3a
template <class Type> 
Real_t MRG32k3a<Type>::fSeed[MRG::vsize] = {12345., 12345., 12345., 12345., 12345., 12345.};

//
// Class Implementation
//   

// Copy constructor
template <typename Type>
VECCORE_ATT_HOST_DEVICE
MRG32k3a<Type>::MRG32k3a(const MRG32k3a<Type> &rng) : VecRNG<MRG32k3a<Type>, Type, MRG32k3a_t<Type> >()
{
  for(int i = 0 ; i < MRG::vsize ; ++i) {
    this->fState->fCg[i] = rng.fState->fCg[i];
    fSeed[i] = rng.fSeed[i]; 
    Bg[i] = rng.Bg[i]; 
  }
}

// Reset stream to the next Stream.
template <typename Type>
VECCORE_ATT_HOST void MRG32k3a<Type>::SetNextStream ()
{
  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) {
    MatVecModM(MRG::A1p127, fSeed, fSeed, MRG::m1);
    MatVecModM(MRG::A2p127, &fSeed[3], &fSeed[3], MRG::m2);
    for (int j = 0; j < MRG::vsize ; ++j) {
      this->fState->fCg[i][j] = Bg[i][j] = fSeed[j];
    }
  }
}

// Scalar specialization of SetNextStream
template <>
VECCORE_ATT_HOST void MRG32k3a<ScalarBackend>::SetNextStream ()
{
  MatVecModM(MRG::A1p127, fSeed, fSeed, MRG::m1);
  MatVecModM(MRG::A2p127, &fSeed[3], &fSeed[3], MRG::m2);

  for (int j = 0; j < MRG::vsize ; ++j) {
    this->fState->fCg[j] = Bg[j] = fSeed[j];
  }
}

// Reset stream to the next substream.
template <typename Type>
VECCORE_ATT_HOST void MRG32k3a<Type>::SetNextSubstream ()
{
  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) {
    MatVecModM(MRG::A1p76, Bg, Bg, MRG::m1);
    MatVecModM(MRG::A2p76, &Bg[3], &Bg[3], MRG::m2);
    for (int j = 0; i < MRG::vsize ; ++i) {
      this->fState->fCg[i][j] = Bg[i][j];
    }
  }
}

// Scalar specialization of SetNextSubstream
template <>
VECCORE_ATT_HOST void MRG32k3a<ScalarBackend>::SetNextSubstream ()
{
  MatVecModM(MRG::A1p76, Bg, Bg, MRG::m1);
  MatVecModM(MRG::A2p76, &Bg[3], &Bg[3], MRG::m2);
  for (int j = 0; j < MRG::vsize ; ++j) {
     this->fState->fCg[j] = Bg[j];
  }
}

// Specialization for the scalar backend to initialize an arrary of states of which size is [blocks*threads]. 
// "states" should be allocated beforehand, but can used for both host and device pointers
template <>
VECCORE_ATT_HOST 
void MRG32k3a<ScalarBackend>::Initialize(MRG32k3a_t<ScalarBackend> *states, int blocks, int threads)
{
  MRG32k3a_t<ScalarBackend>* hstates 
    = (MRG32k3a_t<ScalarBackend> *) malloc (blocks*threads*sizeof(MRG32k3a_t<ScalarBackend>));

  unsigned int nthreads = blocks * threads;
  for (unsigned int tid = 0 ; tid < nthreads ; ++tid) {
    SetNextStream();
    for(size_t j = 0 ; j < MRG::vsize ; ++j) {
      hstates[tid].fCg[j] = this->fState->fCg[j];
    }
  }
#ifdef VECCORE_CUDA
  cudaMemcpy(states, hstates, nthreads*sizeof(MRG32k3a_t<ScalarBackend>), cudaMemcpyHostToDevice);
#else
  memcpy(states, hstates, nthreads*sizeof(MRG32k3a_t<ScalarBackend>));
#endif
  free(hstates);
}

// Print information of the current state
template <typename Type>
VECCORE_ATT_HOST void MRG32k3a<Type>::PrintState() const
{
  for(size_t j = 0 ; j < MRG::vsize ; ++j) {
    std::cout << this->fState->fCg[j] << std::endl;
  }
}

// Set the next seed
template <typename Type>
VECCORE_ATT_HOST_DEVICE void MRG32k3a<Type>::SetSeed(Real_t seed[MRG::vsize]) 
{ 
  for(int i = 0 ; i < MRG::vsize ; ++i) fSeed[i] = seed[i]; 
}

// Kernel to generate the next random number(s) with RetType (based on RngStream::U01d)
template <class Type>
template <class RetType>
inline VECCORE_ATT_HOST_DEVICE typename RetType::Double_v 
MRG32k3a<Type>::Kernel(MRG32k3a_t<Type> *state)
{
  using Double_v = typename RetType::Double_v;

  Double_v k, p1, p2;

  // Component 1 
  p1 = MRG::a12 * state->fCg[1] - MRG::a13n * state->fCg[0];
  k = math::Floor (p1 / MRG::m1);  
  p1 -= k * MRG::m1; 

  Mask_v<Double_v> negative = (p1 < 0.);
  MaskedAssign(p1, negative, p1 + MRG::m1); 

  state->fCg[0] = state->fCg[1];
  state->fCg[1] = state->fCg[2];
  state->fCg[2] = p1;

  // Component 2 
  p2 = MRG::a21 * state->fCg[5] - MRG::a23n * state->fCg[3];
  k = math::Floor (p2 / MRG::m2);
  p2 -= k * MRG::m2;

  negative = (p2 < 0.);
  MaskedAssign(p2, negative, p2 + MRG::m2); 
  
  state->fCg[3] = state->fCg[4];
  state->fCg[4] = state->fCg[5];
  state->fCg[5] = p2;

  // Combination 
  return Blend((p1 > p2),(p1 - p2) * MRG::norm, (p1 - p2 + MRG::m1) * MRG::norm);

  // Extended (53 bits) precision
  // Double_v random =  Blend((p1 > p2),(p1 - p2) * MRG::norm, (p1 - p2 + MRG::m1) * MRG::norm);
  // random *= (1.0 + MRG::fact);
  // return Blend((random < 1.0), random, random - 1.0);
}

// Sepecialization for the scalar method of VectorBackend pRNG class
#ifndef VECPHYS_CUDA
template <>
template <>
inline VECCORE_ATT_HOST_DEVICE ScalarBackend::Double_v 
MRG32k3a<VectorBackend>::Kernel<ScalarBackend>(MRG32k3a_t<VectorBackend> *state)
{
  double k, p1, p2;

  p1 = MRG::a12 * state->fCg[0][1] - MRG::a13n * state->fCg[0][0];
  k = floor (p1 / MRG::m1);
  p1 -= k * MRG::m1;
  if (p1 < 0.0) p1 += MRG::m1;
  state->fCg[0][0] = state->fCg[0][1]; 
  state->fCg[0][1] = state->fCg[0][2]; 
  state->fCg[0][2] = p1;

  // Component 2 
  p2 = MRG::a21 * state->fCg[0][5] - MRG::a23n * state->fCg[0][3];
  k = floor (p2 / MRG::m2);
  p2 -= k * MRG::m2;
  if (p2 < 0.0) p2 += MRG::m2;
  state->fCg[0][3] = state->fCg[0][4]; 
  state->fCg[0][4] = state->fCg[0][5]; 
  state->fCg[0][5] = p2;

  // Combination 
  return ((p1 > p2) ? (p1 - p2) * MRG::norm : (p1 - p2 + MRG::m1) * MRG::norm);
}
#endif

// Sepecialization for the scalar method of ScalarBackend pRNG class
template <>
template <>
inline VECCORE_ATT_HOST_DEVICE ScalarBackend::Double_v 
MRG32k3a<ScalarBackend>::Kernel<ScalarBackend>(MRG32k3a_t<ScalarBackend> *state)
{
  double k, p1, p2;

  p1 = MRG::a12 * state->fCg[1] - MRG::a13n * state->fCg[0];
#if __CUDA_ARCH__ > 0
  k = trunc (fma (p1, MRG::rh1, p1 * MRG::rl1));
  p1 -= k * MRG::m1;
  if (p1 < 0.0) p1 += MRG::m1;
  state->fCg[0] = state->fCg[1]; 
  state->fCg[1] = state->fCg[2]; 
  state->fCg[2] = p1;

  p2 = MRG::a21 * state->fCg[5] - MRG::a23n * state->fCg[3];
  k = trunc (fma (p2, MRG::rh2, p2 * MRG::rl2));
  p2 -= k * MRG::m2;
#else
  k = floor(p1/MRG::m1);
  p1 -= k * MRG::m1;
  if (p1 < 0.0) p1 += MRG::m1;
  state->fCg[0] = state->fCg[1]; 
  state->fCg[1] = state->fCg[2]; 
  state->fCg[2] = p1;

  p2 = MRG::a21 * state->fCg[5] - MRG::a23n * state->fCg[3];
  k = floor(p2/MRG::m2);
  p2 -= k * MRG::m2;
#endif

  if (p2 < 0.0) p2 += MRG::m2;
  state->fCg[3] = state->fCg[4]; state->fCg[4] = state->fCg[5]; state->fCg[5] = p2;

  // Combination 
  return ((p1 > p2) ? (p1 - p2) * MRG::norm : (p1 - p2 + MRG::m1) * MRG::norm);

}

// Return (a*s + c) MOD m; a, s, c and m must be < 2^35
template <class Type>
VECCORE_ATT_HOST double MRG32k3a<Type>::MultModM (double a, double s, double c, double m)
{
  double v;
  long a1;

  v = a * s + c;

  if (v >= MRG::two53 || v <= -MRG::two53) {
    a1 = static_cast<long> (a / MRG::two17);    a -= a1 * MRG::two17;
      v  = a1 * s;
      a1 = static_cast<long> (v / m);     v -= a1 * m;
      v = v * MRG::two17 + a * s + c;
  }

  a1 = static_cast<long> (v / m);
  // in case v < 0)
  if ((v -= a1 * m) < 0.0) return v += m;
  else return v;
}

// Compute the vector v = A*s MOD m. Assume that -m < s[i] < m. Works also when v = s.
template <class Type>
VECCORE_ATT_HOST 
void MRG32k3a<Type>::MatVecModM (const double A[MRG::ndim][MRG::ndim], const double s[MRG::ndim], 
                                 double v[MRG::ndim], double m)
{
  int i;
  double x[MRG::ndim];  

  for (i = 0; i < MRG::ndim ; ++i) {
    x[i] = MultModM (A[i][0], s[0], 0.0, m);
    x[i] = MultModM (A[i][1], s[1], x[i], m);
    x[i] = MultModM (A[i][2], s[2], x[i], m);
  }
  for (i = 0; i < MRG::ndim ; ++i) v[i] = x[i];
}

} // end namespace impl
} // end namespace vecphys

#endif
