#ifndef GEANTV_COMMONBENCHMARK_H
#define GEANTV_COMMONBENCHMARK_H

#include <benchmark/benchmark.h>
#include "Common.h"

#include <functional>

void ScalarModelBenchmark(benchmark::State &state, EMModel *model,
                          std::function<void(std::vector<LightTrack> &, int)> preparePrimaries, TaskData *td,
                          int basketSize)
{
  std::vector<LightTrack> primaries;

  for (auto _ : state) {
    state.PauseTiming();

    preparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      model->SampleSecondaries(primaries[i], td);
    }
  }

  benchmark::DoNotOptimize(primaries.data());
}

void VectorModelBenchmark(benchmark::State &state, EMModel *model,
                          std::function<void(LightTrack_v &, int)> preparePrimaries, TaskData *td, int basketSize)
{
  LightTrack_v primaries;

  for (auto _ : state) {
    state.PauseTiming();

    preparePrimaries(primaries, basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    model->SampleSecondaries(primaries, td);
  }

  benchmark::DoNotOptimize(&primaries);
}

#define EMMODEL_REGISTER_SCALAR_BENCHMARK(Name, Model, PrimaryGenFunction, TaskData, BasketSize) \
  benchmark::RegisterBenchmark((Name), [&](benchmark::State &st) {                               \
    ScalarModelBenchmark(st, (Model), (PrimaryGenFunction), (TaskData), (BasketSize));           \
  });

#define EMMODEL_REGISTER_VECTOR_BENCHMARK(Name, Model, PrimaryGenFunction, TaskData, BasketSize) \
  benchmark::RegisterBenchmark((Name), [&](benchmark::State &st) {                               \
    VectorModelBenchmark(st, (Model), (PrimaryGenFunction), (TaskData), (BasketSize));           \
  });

#endif // GEANTV_COMMONBENCHMARK_H
