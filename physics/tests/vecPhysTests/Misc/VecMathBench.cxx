#include <benchmark/benchmark.h>
#include <Geant/math_wrappers.h>

#include "Geant/RngWrapper.h"
#include "Geant/VectorTypes.h"
using geant::Double_v;
using geant::kVecAlignD;
using geant::kVecLenD;

const int kN = 256;

static void ScalarExp(benchmark::State &state)
{
  geant::RngWrapper wrapper;
  std::vector<double> rnd;
  for (int i = 0; i < state.range(0); ++i) {
    rnd.push_back(wrapper.uniform());
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); ++i) {
      sum += Math::Exp(rnd[i]);
    }
    benchmark::DoNotOptimize(sum);
  }
}
BENCHMARK(ScalarExp)->Arg(kN);

static void VectorExp(benchmark::State &state)
{
  geant::RngWrapper wrapper;

  double *rnd = (double *)vecCore::AlignedAlloc(kVecAlignD, state.range(0) * sizeof(double));
  for (int i = 0; i < state.range(0); ++i) {
    rnd[i] = wrapper.uniform();
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); i += kVecLenD) {
      Double_v tmp;
      vecCore::Load(tmp, &rnd[i]);
      tmp = Math::Exp(tmp);
      sum += vecCore::ReduceAdd(tmp);
    }
    benchmark::DoNotOptimize(sum);
  }
  vecCore::AlignedFree(rnd);
}
BENCHMARK(VectorExp)->Arg(kN);

static void ScalarLog(benchmark::State &state)
{
  geant::RngWrapper wrapper;
  std::vector<double> rnd;
  for (int i = 0; i < state.range(0); ++i) {
    rnd.push_back(wrapper.uniform());
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); ++i) {
      sum += Math::Log(rnd[i]);
    }
    benchmark::DoNotOptimize(sum);
  }
}
BENCHMARK(ScalarLog)->Arg(kN);

static void VectorLog(benchmark::State &state)
{
  geant::RngWrapper wrapper;

  double *rnd = (double *)vecCore::AlignedAlloc(kVecAlignD, state.range(0) * sizeof(double));
  for (int i = 0; i < state.range(0); ++i) {
    rnd[i] = wrapper.uniform();
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); i += kVecLenD) {
      Double_v tmp;
      vecCore::Load(tmp, &rnd[i]);
      tmp = Math::Log(tmp);
      sum += vecCore::ReduceAdd(tmp);
    }
    benchmark::DoNotOptimize(sum);
  }
  vecCore::AlignedFree(rnd);
}
BENCHMARK(VectorLog)->Arg(kN);

static void ScalarSqrt(benchmark::State &state)
{
  geant::RngWrapper wrapper;
  std::vector<double> rnd;
  for (int i = 0; i < state.range(0); ++i) {
    rnd.push_back(wrapper.uniform());
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); ++i) {
      sum += std::sqrt(rnd[i]);
    }
    benchmark::DoNotOptimize(sum);
  }
}
BENCHMARK(ScalarSqrt)->Arg(kN);

static void VectorSqrt(benchmark::State &state)
{
  geant::RngWrapper wrapper;

  double *rnd = (double *)vecCore::AlignedAlloc(kVecAlignD, state.range(0) * sizeof(double));
  for (int i = 0; i < state.range(0); ++i) {
    rnd[i] = wrapper.uniform();
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); i += kVecLenD) {
      Double_v tmp;
      vecCore::Load(tmp, &rnd[i]);
      tmp = Math::Sqrt(tmp);
      sum += vecCore::ReduceAdd(tmp);
    }
    benchmark::DoNotOptimize(sum);
  }
  vecCore::AlignedFree(rnd);
}
BENCHMARK(VectorSqrt)->Arg(kN);

static void ScalarDiv(benchmark::State &state)
{
  geant::RngWrapper wrapper;
  std::vector<double> rnd;
  for (int i = 0; i < state.range(0); ++i) {
    rnd.push_back(wrapper.uniform());
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); ++i) {
      double tmp = rnd[i];
      tmp        = 1.0 / tmp;
      tmp        = 1.0 / tmp;
      tmp        = 1.0 / tmp;
      tmp        = 1.0 / tmp;
      sum += tmp;
    }
    benchmark::DoNotOptimize(sum);
  }
}
BENCHMARK(ScalarDiv)->Arg(kN);

static void VectorDiv(benchmark::State &state)
{
  geant::RngWrapper wrapper;

  double *rnd = (double *)vecCore::AlignedAlloc(kVecAlignD, state.range(0) * sizeof(double));
  for (int i = 0; i < state.range(0); ++i) {
    rnd[i] = wrapper.uniform();
  }
  for (auto _ : state) {
    Double_v sum = 0.0;
    for (int i = 0; i < state.range(0); i += kVecLenD) {
      Double_v tmp;
      vecCore::Load(tmp, &rnd[i]);
      tmp = 1.0 / tmp;
      tmp = 1.0 / tmp;
      tmp = 1.0 / tmp;
      tmp = 1.0 / tmp;
      sum += tmp;
    }
    benchmark::DoNotOptimize(sum);
  }
  vecCore::AlignedFree(rnd);
}
BENCHMARK(VectorDiv)->Arg(kN);

static void ScalarSin(benchmark::State &state)
{
  geant::RngWrapper wrapper;
  std::vector<double> rnd;
  for (int i = 0; i < state.range(0); ++i) {
    rnd.push_back(wrapper.uniform());
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); ++i) {
      sum += Math::Sin(rnd[i]);
    }
    benchmark::DoNotOptimize(sum);
  }
}
BENCHMARK(ScalarSin)->Arg(kN);

static void VectorSin(benchmark::State &state)
{
  geant::RngWrapper wrapper;

  double *rnd = (double *)vecCore::AlignedAlloc(kVecAlignD, state.range(0) * sizeof(double));
  for (int i = 0; i < state.range(0); ++i) {
    rnd[i] = wrapper.uniform();
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); i += kVecLenD) {
      Double_v tmp;
      vecCore::Load(tmp, &rnd[i]);
      tmp = Math::Sin(tmp);
      sum += vecCore::ReduceAdd(tmp);
    }
    benchmark::DoNotOptimize(sum);
  }
  vecCore::AlignedFree(rnd);
}
BENCHMARK(VectorSin)->Arg(kN);

static void ScalarCos(benchmark::State &state)
{
  geant::RngWrapper wrapper;
  std::vector<double> rnd;
  for (int i = 0; i < state.range(0); ++i) {
    rnd.push_back(wrapper.uniform());
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); ++i) {
      sum += Math::Cos(rnd[i]);
    }
    benchmark::DoNotOptimize(sum);
  }
}
BENCHMARK(ScalarCos)->Arg(kN);

static void VectorCos(benchmark::State &state)
{
  geant::RngWrapper wrapper;

  double *rnd = (double *)vecCore::AlignedAlloc(kVecAlignD, state.range(0) * sizeof(double));
  for (int i = 0; i < state.range(0); ++i) {
    rnd[i] = wrapper.uniform();
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); i += kVecLenD) {
      Double_v tmp;
      vecCore::Load(tmp, &rnd[i]);
      tmp = Math::Cos(tmp);
      sum += vecCore::ReduceAdd(tmp);
    }
    benchmark::DoNotOptimize(sum);
  }
  vecCore::AlignedFree(rnd);
}
BENCHMARK(VectorCos)->Arg(kN);

static void ScalarSinCos(benchmark::State &state)
{
  geant::RngWrapper wrapper;
  std::vector<double> rnd;
  for (int i = 0; i < state.range(0); ++i) {
    rnd.push_back(wrapper.uniform());
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); ++i) {
      sum += Math::Sin(rnd[i]);
      sum += Math::Cos(rnd[i]);
    }
    benchmark::DoNotOptimize(sum);
  }
}
BENCHMARK(ScalarSinCos)->Arg(kN);

static void VectorSinCos(benchmark::State &state)
{
  geant::RngWrapper wrapper;

  double *rnd = (double *)vecCore::AlignedAlloc(kVecAlignD, state.range(0) * sizeof(double));
  for (int i = 0; i < state.range(0); ++i) {
    rnd[i] = wrapper.uniform();
  }
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); i += kVecLenD) {
      Double_v tmp;
      vecCore::Load(tmp, &rnd[i]);
      Double_v s, c;
      Math::SinCos(tmp, s, c);
      sum += vecCore::ReduceAdd(s);
      sum += vecCore::ReduceAdd(c);
    }
    benchmark::DoNotOptimize(sum);
  }
  vecCore::AlignedFree(rnd);
}
BENCHMARK(VectorSinCos)->Arg(kN);

BENCHMARK_MAIN();
