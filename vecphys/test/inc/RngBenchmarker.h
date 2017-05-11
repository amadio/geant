#ifndef RngBenchmarker_H
#define RngBenchmarker_H 1

#include "base/VecPhys.h"

namespace vecphys {

class RngBenchmarker {

public:
  RngBenchmarker();
  ~RngBenchmarker();

  int RunBenchmark();

  void SetNSample(const int nsample) { fNSample = nsample; }
  void SetRepetition(const unsigned repetition) { fRepetition = repetition; }

private:
  int RunBenchmarkRng();

  void RunTest();
  void RunScalar();
  void RunVector();
  void RunNState();

#ifdef VECPHYS_MKL
  void RunMKLVSL();
#endif

#ifdef VECPHYS_CUDA
  void RunCuda();
#endif

private:
  int fNSample;
  unsigned fRepetition;
  int fVerbosity;

};

} // end namespace vecphys

#endif
