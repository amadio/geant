#include "RngBenchmarker.h"
#include "RngBenchmarker_cpu.h"

#include "base/Stopwatch.h"

using vecgeom::Stopwatch;

namespace vecphys {

enum RngIndex { kNullRng = -1, kMRG32k3a, kThreefry, kPhilox, kNumberRng };
static const char *RngName[kNumberRng] = {"MRG32k3a", "Threefry", "Philox  "};

#ifdef VECPHYS_MKL
  enum VslIndex { kNullVsl = -1, kMRG32K3A, kMCG59, kMT17739, kSFMT17739, kGFSR250, kSOBOL32, kNumberVsl };
  static const char *VslName[kNumberVsl] = {"MRG32K3A", "MCG59   ", "MT19937 ", "SMT19937","GFSR250 ","SOBOL32 "};
#endif

RngBenchmarker::RngBenchmarker()
  : fNSample(100), fRepetition(1)
{
}

RngBenchmarker::~RngBenchmarker()
{
}

int RngBenchmarker::RunBenchmark()
{
  printf(" Run RngBenchmarker with arguments: %d %d\n", fNSample, fRepetition);
  printf(" NSample              = %d\n", fNSample);
  printf(" NRepetitions         = %d\n", fRepetition);
  printf(" VectorSize<Double_v> = %lu\n", VectorSize<Double_v>());

  int errorcode = 0;
  errorcode += RunBenchmarkRng();
  return (errorcode) ? 1 : 0;
}

int RngBenchmarker::RunBenchmarkRng()
{
  if (fVerbosity > 0) {
  }

  //  int mismatches = 0;

  RunTest();
  RunScalar();
  RunVector();
  RunNState();

#ifdef VECPHYS_MKL
  RunMKLVSL();
#endif

#ifdef VECPHYS_CUDA
    RunCuda();
#endif

    //  return mismatches;
  return 0;
}

void RngBenchmarker::RunTest()
{
  // test GNU rand() as a reference

  Real_t elapsedTotal = 0.;
  double resultTotal = 0.;
  static Stopwatch timer;

  for (unsigned r = 0; r < fRepetition; ++r) {

    Real_t elapsedTime = 0.;
    double result = 0.;

    timer.Start();

    for (int i = 0; i < fNSample ; ++i) {
      result += (double)rand()/RAND_MAX;
    }

    elapsedTime = timer.Stop();

    elapsedTotal += elapsedTime;
    resultTotal += result;
  }

  printf(" %s std::rand()   Total time = %6.3f sec CheckSum = %g\n", 
	 "TestRand ", elapsedTotal, resultTotal);
}

void RngBenchmarker::RunScalar()
{
  Real_t elapsedTotal[kNumberRng];
  Real_t elapsedT;

  double resultTotal[kNumberRng];
  double result;

  for (int k = 0; k < kNumberRng ; ++k) {
    elapsedTotal[k] = 0.;
    resultTotal[k] = 0.;
  }
  for (unsigned r = 0; r < fRepetition; ++r) {
    for (int k = 0; k < kNumberRng ; ++k) {
      elapsedT = 0.0;
      result = 0.0;
      elapsedT = ScalarKernelFunc[k](fNSample,result);
      elapsedTotal[k] += elapsedT;
      resultTotal[k] += result;
    }
  }

  for (int k = 0; k < kNumberRng; ++k) {
    printf(" %s  ScalarBackend Total time = %6.3f sec CheckSum = %g\n", 
	   RngName[k], elapsedTotal[k], resultTotal[k]);
  }
}

void RngBenchmarker::RunVector()
{
  Real_t elapsedTotal[kNumberRng];
  Real_t elapsedT;

  double resultTotal[kNumberRng];
  double result;

  for (int k = 0; k < kNumberRng ; ++k) {
    elapsedTotal[k] = 0.;
    resultTotal[k] = 0.;
  }
  for (unsigned r = 0; r < fRepetition; ++r) {
    for (int k = 0; k < kNumberRng ; ++k) {
      elapsedT = 0.0;
      result = 0.0;
      elapsedT = VectorKernelFunc[k](fNSample,result);
      elapsedTotal[k] += elapsedT;
      resultTotal[k] += result;
    }
  }

  for (int k = 0; k < kNumberRng; ++k) {
    printf(" %s  VectorBackend Total time = %6.3f sec CheckSum = %g\n", 
	   RngName[k], elapsedTotal[k], resultTotal[k]);
  }

}

void RngBenchmarker::RunNState()
{
  Real_t elapsedTotal[kNumberRng];
  Real_t elapsedT;

  double resultTotal[kNumberRng];
  double result;

  for (int k = 0; k < kNumberRng ; ++k) {
    elapsedTotal[k] = 0.;
    resultTotal[k] = 0.;
  }
  for (unsigned r = 0; r < fRepetition; ++r) {
    for (int k = 0; k < kNumberRng ; ++k) {
      elapsedT = 0.0;
      result = 0.0;
      elapsedT = StateKernelFunc[k](fNSample,result);
      elapsedTotal[k] += elapsedT;
      resultTotal[k] += result;
    }
  }

  for (int k = 0; k < kNumberRng; ++k) {
    printf(" %s  ScalarNstates Total time = %6.3f sec CheckSum = %g\n", 
	   RngName[k], elapsedTotal[k], resultTotal[k]);
  }
}

#ifdef VECPHYS_MKL
void RngBenchmarker::RunMKLVSL()
{
  Real_t elapsedTotal[kNumberVsl];
  Real_t elapsedT;

  double resultTotal[kNumberVsl];
  double result;

  for (int k = 0; k < kNumberVsl ; ++k) {
    elapsedTotal[k] = 0.;
    resultTotal[k] = 0.;
  }
  for (unsigned r = 0; r < fRepetition; ++r) {
    for (unsigned int k = 0; k < kNumberVsl ; ++k) {
      elapsedT = 0.0;
      result = 0.0;
      elapsedT = VSLKernelFunc[k](fNSample,result);
      elapsedTotal[k] += elapsedT;
      resultTotal[k] += result;
    }
  }

  for (int k = 0; k < kNumberVsl; ++k) {
    printf(" %s  Intel-MKL/VSL Total time = %6.3f sec CheckSum = %g\n", 
	   VslName[k], elapsedTotal[k], resultTotal[k]);
  }
}
#endif

} // end namespace vecphys
