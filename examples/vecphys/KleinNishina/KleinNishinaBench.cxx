#include <benchmark/benchmark.h>
#include "KleinNishinaTestCommon.h"

#include "Geant/VecRngWrapper.h"

static void KleinNishinaScalarAliasBaseline(benchmark::State &state)
{
  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(true);

  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    PreparePrimaries(primaries, basketSize);
    Td->fPhysicsData->ClearSecondaries();
  }

  benchmark::DoNotOptimize(primaries.data());

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaScalarAliasBaseline)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void KleinNishinaScalarAlias(benchmark::State &state)
{
  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(true);

  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    PreparePrimaries(primaries, basketSize);
    Td->fPhysicsData->ClearSecondaries();
    for (int i = 0; i < basketSize; ++i) {
      kNish->SampleSecondaries(primaries[i], Td);
    }
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void KleinNishinaVectorAliasBaseline(benchmark::State &state)
{
  VecKleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(true);

  auto Td = PrepareTaskData();

  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);

    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
  }
  benchmark::DoNotOptimize(&primaries);

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaVectorAliasBaseline)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void KleinNishinaVectorAlias(benchmark::State &state)
{
  VecKleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(true);

  auto Td = PrepareTaskData();

  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);

    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
    kNish->SampleSecondariesVector(primaries, Td);
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaVectorAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void KleinNishinaScalarRej(benchmark::State &state)
{

  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(false);

  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    state.ResumeTiming();
    Td->fPhysicsData->ClearSecondaries();
    for (int i = 0; i < basketSize; ++i) {
      kNish->SampleSecondaries(primaries[i], Td);
    }
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

BENCHMARK_MAIN();
