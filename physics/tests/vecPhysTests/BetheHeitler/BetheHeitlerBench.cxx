#include <benchmark/benchmark.h>
#include "BetheHeitlerTestCommon.h"

#include "Geant/VecRngWrapper.h"

BetheHeitlerPairModel *rpModelRej;
VecBetheHeitlerPairModel *rpModelRejVec;
BetheHeitlerPairModel *rpModelAlias;
VecBetheHeitlerPairModel *rpModelAliasVec;
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
      rpModelAlias->SampleSecondaries(primaries[i], td);
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

    rpModelAliasVec->SampleSecondariesVector(primaries, td);
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
      rpModelRej->SampleSecondaries(primaries[i], td);
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

    rpModelRejVec->SampleSecondariesVector(primaries, td);
  }
}
BENCHMARK(BetheHeitlerVectorRej)->RangeMultiplier(2)->Range(kMinBasket, kMaxBasket);

int main(int argc, char **argv)
{
  PrepareWorld();
  rpModelAlias    = PrepareBHModel(true);
  rpModelRej      = PrepareBHModel(false);
  rpModelAliasVec = PrepareVecBHModel(true);
  rpModelRejVec   = PrepareVecBHModel(false);
  td              = PrepareTaskData();

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
  ::benchmark::RunSpecifiedBenchmarks();

  CleanTaskData(td);
}
