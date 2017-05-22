#include "base/Stopwatch.h"
using vecgeom::Stopwatch;

#include "base/VecPhys.h"

#include "rng/MRG32k3a.h"
#include "rng/Threefry.h"
#include "rng/Philox.h"

#ifdef VECPHYS_MKL
#include "mkl_vsl.h"
#include <omp.h>
#endif

namespace vecphys {

// Scalar

double ScalarMRG32k3a(int nsample, double& result)
{
  // Scalar MRG32k3a
  static vecphys::cxx::MRG32k3a<ScalarBackend> rng;
  rng.Initialize();

  static Stopwatch timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample ; ++i) {
    sum += rng.Uniform<ScalarBackend>();
  }

  elapsedTime = timer.Stop();
  result = sum;

  return elapsedTime;
}

double ScalarThreefry(int nsample, double& result)
{
  // Scalar Threefry
  static vecphys::cxx::Threefry<ScalarBackend> rng;
  rng.Initialize();

  static Stopwatch timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample ; ++i) {
    sum += rng.Uniform<ScalarBackend>();
  }

  elapsedTime = timer.Stop();
  result = sum;

  return elapsedTime;
}

double ScalarPhilox(int nsample, double& result)
{
  // Scalar Philox - use scalar method of VectorBackend (temporarily)

  static vecphys::cxx::Philox<ScalarBackend> rng;
  rng.Initialize();

  static Stopwatch timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample ; ++i) {
    sum += rng.Uniform<ScalarBackend>();
  }

  elapsedTime = timer.Stop();
  result = sum;

  vecphys::cxx::Philox<ScalarBackend> rng1 = rng;

  return elapsedTime;
}

// Vector

double VectorMRG32k3a(int nsample, double& result)
{
  // Scalar MRG32k3a
  using Double_v = typename VectorBackend::Double_v;

  static vecphys::cxx::MRG32k3a<VectorBackend> rng;
  rng.Initialize();

  static Stopwatch timer;
  double elapsedTime = 0.;

  Double_v sum = 0;

  timer.Start();

  for (int i = 0; i < nsample/2 ; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Stop();
  for (unsigned int i = 0; i < VectorSize<Double_v>() ; ++i) result += sum[i];

  return elapsedTime;
}

double VectorThreefry(int nsample, double& result)
{
  // Vector Threefry
  using Double_v = typename VectorBackend::Double_v;

  static vecphys::cxx::Threefry<VectorBackend> rng;
  rng.Initialize();

  static Stopwatch timer;
  double elapsedTime = 0.;

  Double_v sum = 0;

  timer.Start();

  for (int i = 0; i < nsample/2 ; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Stop();
  for (unsigned int i = 0; i < VectorSize<Double_v>() ; ++i) result += sum[i];

  return elapsedTime;
}

double VectorPhilox(int nsample, double& result)
{
  // Vector Philox
  using Double_v = typename VectorBackend::Double_v;

  static vecphys::cxx::Philox<VectorBackend> rng;
  rng.Initialize();

  static Stopwatch timer;
  double elapsedTime = 0.;

  Double_v sum = 0;

  timer.Start();

  for (int i = 0; i < nsample/2 ; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Stop();
  for (unsigned int i = 0; i < VectorSize<Double_v>() ; ++i) result += sum[i];

  return elapsedTime;
}

double StateMRG32k3a(int nsample, double& result)
{
  // Scalar MRG32k3a
  static vecphys::cxx::MRG32k3a<ScalarBackend> rng;

  int theNBlocks = 26;
  int theNThreads = 192;

  MRG32k3a_t<ScalarBackend>* hstates 
    = (MRG32k3a_t<ScalarBackend> *) malloc (theNBlocks*theNThreads*sizeof(MRG32k3a_t<ScalarBackend>));
  rng.Initialize(hstates,theNBlocks,theNThreads);

  static Stopwatch timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for(int i = 0 ; i < theNBlocks*theNThreads ; ++i) {
    for (int j = 0; j < nsample/(theNBlocks*theNThreads) ; ++j) {
      sum += rng.Uniform<ScalarBackend>(&hstates[i]);
    }
  }

  elapsedTime = timer.Stop();
  result = sum;

  return elapsedTime;
}

double StateThreefry(int nsample, double& result) 
{
  // Scalar Threefry
  static vecphys::cxx::Threefry<ScalarBackend> rng;

  int theNBlocks = 26;
  int theNThreads = 192;

  Threefry_t<ScalarBackend>* hstates 
    = (Threefry_t<ScalarBackend> *) malloc (theNBlocks*theNThreads*sizeof(Threefry_t<ScalarBackend>));
  rng.Initialize(hstates,theNBlocks,theNThreads);

  static Stopwatch timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for(int i = 0 ; i < theNBlocks*theNThreads ; ++i) {
    for (int j = 0; j < nsample/(theNBlocks*theNThreads) ; ++j) {
      sum += rng.Uniform<ScalarBackend>(&hstates[i]);
    }
  }

  elapsedTime = timer.Stop();
  result = sum;

  return elapsedTime;
}

double StatePhilox(int nsample, double& result) 
{
  // Scalar Philox
  //  static vecphys::cxx::VecRNG<Philox<ScalarBackend>, ScalarBackend, Philox_t<ScalarBackend>> rng;
  static vecphys::cxx::Philox<ScalarBackend> rng;

  int theNBlocks = 26;
  int theNThreads = 192;

  Philox_t<ScalarBackend>* hstates
    = (Philox_t<ScalarBackend> *) malloc (theNBlocks*theNThreads*sizeof(Philox_t<ScalarBackend>));
  rng.Initialize(hstates,theNBlocks,theNThreads);

  static Stopwatch timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for(int i = 0 ; i < theNBlocks*theNThreads ; ++i) {
    for (int j = 0; j < nsample/(theNBlocks*theNThreads) ; ++j) {
      sum += rng.Uniform<ScalarBackend>(&hstates[i]);
    }
  }

  elapsedTime = timer.Stop();
  result = sum;

  return elapsedTime;
}

// Intel VMK/VSL
#ifdef VECPHYS_MKL

double VSLRngTest(int nsample, double& result, int rngType)
{
  static Stopwatch timer;
  double elapsedTime = 0.;

  const int N = 1;

  double r[N];
  double a=0.0;
  double b=1.0;

  double sum = 0;

  // Initialize
  unsigned int seed = 7777777;
  VSLStreamStatePtr stream;
  vslNewStream(&stream, rngType, (MKL_INT)seed);

  timer.Start();

  // Call RNG
  for (int i = 0; i < nsample/N ; ++i) {
    vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, N, r, a, b);
    for(int j = 0 ; j < N ; ++j) sum += r[j];
  }
  elapsedTime = timer.Stop();

  result = sum;

  // Deinitialize
  vslDeleteStream( &stream );

  return elapsedTime;
}

#endif

} // end namespace vecphys
