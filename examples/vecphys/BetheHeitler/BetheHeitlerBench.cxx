#include <benchmark/benchmark.h>
#include "BetheHeitlerTestCommon.h"

#include "Geant/VecRngWrapper.h"

BetheHeitlerPairModel *bHModelRej;
VecBetheHeitlerPairModel *bHVecModelRej;
BetheHeitlerPairModel *bHModelAlias;
VecBetheHeitlerPairModel *bHVecModelAlias;
TaskData *td;

static void BetheHeitlerScalarAlias(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      bHModelAlias->SampleSecondaries(primaries[i], td);
    }
  }
}
BENCHMARK(BetheHeitlerScalarAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void BetheHeitlerVectorAlias(benchmark::State &state)
{
  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    bHVecModelAlias->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(BetheHeitlerVectorAlias)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void BetheHeitlerScalarRej(benchmark::State &state)
{
  std::vector<LightTrack> primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    td->fPhysicsData->ClearSecondaries();

    state.ResumeTiming();

    for (int i = 0; i < basketSize; ++i) {
      bHModelRej->SampleSecondaries(primaries[i], td);
    }
  }
}
BENCHMARK(BetheHeitlerScalarRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

static void BetheHeitlerVectorRej(benchmark::State &state)
{
  LightTrack_v primaries;

  int basketSize = state.range(0);

  for (auto _ : state) {
    state.PauseTiming();

    PreparePrimaries(primaries, basketSize);
    primaries.SetNtracks(basketSize);
    td->fPhysicsData->GetSecondarySOA().ClearTracks();

    state.ResumeTiming();

    bHVecModelRej->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(BetheHeitlerVectorRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

int main(int argc, char **argv)
{
  PrepareWorld();
  bHModelAlias    = PrepareBHModel(true);
  bHModelRej      = PrepareBHModel(false);
  bHVecModelAlias = PrepareVecBHModel(true);
  bHVecModelRej   = PrepareVecBHModel(false);
  td              = PrepareTaskData();

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td);
}
