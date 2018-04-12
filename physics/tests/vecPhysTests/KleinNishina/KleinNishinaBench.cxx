#include <benchmark/benchmark.h>
#include "KleinNishinaTestCommon.h"

#include "Geant/VecRngWrapper.h"

static void KleinNishinaScalarAlias(benchmark::State &state)
{
  KleinNishinaComptonModel *kNish = PrepareKnishinaModel(true);

  auto Td = PrepareTaskData();
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    Td->fPhysicsData->ClearSecondaries();
    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      kNish->SampleSecondaries(primaries[i], Td);
    }
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void KleinNishinaVectorAlias(benchmark::State &state)
{
  VecKleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(true);

  auto Td = PrepareTaskData();

  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    state.ResumeTiming();

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
    Td->fPhysicsData->ClearSecondaries();
    state.ResumeTiming();
    for (int i = 0; i < basketSize; ++i) {
      kNish->SampleSecondaries(primaries[i], Td);
    }
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void KleinNishinaVectorRej(benchmark::State &state)
{
  VecKleinNishinaComptonModel *kNish = PrepareVecKnishinaModel(false);

  auto Td = PrepareTaskData();

  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();
    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    state.ResumeTiming();

    Td->fPhysicsData->GetSecondarySOA().ClearTracks();
    kNish->SampleSecondariesVector(primaries, Td);
  }

  delete kNish;
  CleanTaskData(Td);
}
BENCHMARK(KleinNishinaVectorRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

BENCHMARK_MAIN();
