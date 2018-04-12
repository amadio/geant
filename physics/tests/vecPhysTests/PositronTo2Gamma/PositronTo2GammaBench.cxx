#include <benchmark/benchmark.h>
#include "PositronTo2GammaTestCommon.h"

PositronTo2GammaModel *baseAlias;
PositronTo2GammaModel *baseRej;
VecPositronTo2GammaModel *vectorAlias;
VecPositronTo2GammaModel *vectorRej;

static void PosTo2GammaScalarAlias(benchmark::State &state)
{

  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    Td->fPhysicsData->ClearSecondaries();
    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      baseAlias->SampleSecondaries(primaries[i], Td);
    }
  }

  CleanTaskData(Td);
}
BENCHMARK(PosTo2GammaScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void PosTo2GammaVectorAlias(benchmark::State &state)
{
  auto Td = PrepareTaskData();

  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    state.ResumeTiming();

    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
    vectorAlias->SampleSecondariesVector(primaries, Td);
  }

  CleanTaskData(Td);
}
BENCHMARK(PosTo2GammaVectorAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void PosTo2GammaScalarRej(benchmark::State &state)
{
  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    Td->fPhysicsData->ClearSecondaries();
    state.ResumeTiming();
    for (int i = 0; i < basketSize; ++i) {
      baseRej->SampleSecondaries(primaries[i], Td);
    }
  }

  CleanTaskData(Td);
}
BENCHMARK(PosTo2GammaScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void PosTo2GammaVectorRej(benchmark::State &state)
{
  auto Td = PrepareTaskData();

  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    state.ResumeTiming();

    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
    vectorRej->SampleSecondariesVector(primaries, Td);
  }

  CleanTaskData(Td);
}
BENCHMARK(PosTo2GammaVectorRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

int main(int argc, char **argv)
{
  baseAlias   = PrepareAnihilModel(true);
  baseRej     = PrepareAnihilModel(false);
  vectorAlias = PrepareVecAnihilModel(true);
  vectorRej   = PrepareVecAnihilModel(false);

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();
}
