#include <benchmark/benchmark.h>

#include "Geant/RngWrapper.h"
#include "Geant/VectorPhysicsTypes.h"

static void ScalarRNG(benchmark::State &state)
{
  geant::RngWrapper wrapper;
  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < state.range(0); ++i) {
      sum += wrapper.uniform();
    }
    benchmark::DoNotOptimize(sum);
  }
}
BENCHMARK(ScalarRNG)->Range(1 << 3, 1 << 9)->RangeMultiplier(2);

static void VectorRNG(benchmark::State &state)
{
  geant::RngWrapper wrapper;
  for (auto _ : state) {
    PhysDV sumV = 0.0;
    for (int i = 0; i < state.range(0) / kPhysDVWidth; ++i) {
      sumV += wrapper.uniformV();
    }
    double sum = vecCore::ReduceAdd(sum);
    benchmark::DoNotOptimize(sum);
  }
}
BENCHMARK(VectorRNG)->Range(1 << 3, 1 << 9)->RangeMultiplier(2);

BENCHMARK_MAIN();
